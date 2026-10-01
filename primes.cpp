// primes.cpp - Find the first N primes (default: one billion) with a
// multithreaded, segmented Sieve of Eratosthenes using mod-30 wheel
// factorization.
//
// Design
// ------
//  * Mod-30 wheel: only numbers coprime to 2, 3 and 5 are stored. Each byte
//    covers 30 consecutive integers and its 8 bits stand for the residues
//    {1, 7, 11, 13, 17, 19, 23, 29}. That is 8 bits per 30 numbers, i.e. the
//    sieve touches 73% less memory than a plain bit sieve and never looks at
//    a multiple of 2, 3 or 5.
//  * Wheel-driven crossing off: a sieving prime p only strikes p*q where q is
//    itself coprime to 30. The byte step and bit mask for every (p mod 30,
//    q mod 30) pair are precomputed, and the 8 strikes of one wheel turn are
//    unrolled (one turn advances exactly p bytes).
//  * Pre-sieving: multiples of 7, 11, 13 and 17 are removed by copying a
//    repeating 17017-byte pattern instead of being crossed off.
//  * Segmentation: the range is processed in cache-sized segments so the hot
//    working set stays in L1/L2.
//  * Threads grab segments dynamically; each segment is independent, so the
//    only shared state is an atomic counter and a per-segment count array.
//  * The N-th prime is located from per-segment prime counts, after which
//    only the one segment that contains it is sieved again and scanned.
//
// Build:  g++ -O3 -march=native -std=c++17 -pthread primes.cpp -o primes
// Usage:  ./primes [N] [--threads T] [--segment KB] [--print FILE]
//
//   N            how many primes to find (default 1000000000)
//   --threads    worker threads (default: hardware concurrency)
//   --segment    segment size in KiB (default 256)
//   --print      also write all N primes, one per line, to FILE ("-" = stdout)

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

namespace {

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kResidues[8] = {1, 7, 11, 13, 17, 19, 23, 29};
// Gap from residue k to residue k+1 (the last one wraps to 31 = 30 + 1).
constexpr u32 kGaps[8] = {6, 4, 2, 4, 2, 4, 6, 2};

// Pre-sieve primes and the pattern period in bytes (7*11*13*17).
constexpr u32 kPreSievePrimes[4] = {7, 11, 13, 17};
constexpr u32 kPatternBytes = 7 * 11 * 13 * 17;
constexpr u32 kFirstSievingPrime = 19;

struct Tables {
  int residue_index[30];  // residue -> bit index, or -1 if not coprime to 30
  u8 mask[8][8];          // [p idx][q idx] -> AND mask that clears p*q's bit
  u32 carry[8][8];        // [p idx][q idx] -> floor(rp*rq'/30) - floor(rp*rq/30)
  std::vector<u8> pattern;

  Tables() {
    std::fill(std::begin(residue_index), std::end(residue_index), -1);
    for (int k = 0; k < 8; ++k) residue_index[kResidues[k]] = k;

    // With p = 30a + rp and q = 30b + rq, the byte holding p*q is
    //   30ab + a*rq + b*rp + floor(rp*rq/30)
    // so stepping q to the next wheel residue moves the byte by
    //   a*gap + floor(rp*rq_next/30) - floor(rp*rq/30).
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        u32 rp = kResidues[i], rq = kResidues[j];
        u32 rq_next = rq + kGaps[j];
        mask[i][j] = static_cast<u8>(~(1u << residue_index[(rp * rq) % 30]));
        carry[i][j] = (rp * rq_next) / 30 - (rp * rq) / 30;
      }
    }

    // Pattern with every multiple of 7, 11, 13, 17 (including the primes
    // themselves) removed; it repeats every 30*17017 integers.
    pattern.assign(kPatternBytes, 0xFF);
    for (u32 p : kPreSievePrimes) {
      for (u64 m = p; m < 30ull * kPatternBytes; m += 2 * p) {
        int r = residue_index[m % 30];
        if (r >= 0) pattern[m / 30] &= static_cast<u8>(~(1u << r));
      }
    }
  }
};

const Tables& tables() {
  static const Tables t;
  return t;
}

// Primes in [kFirstSievingPrime, limit] via a simple sieve.
std::vector<u32> sieving_primes(u32 limit) {
  std::vector<u8> composite(limit + 1, 0);
  std::vector<u32> out;
  for (u32 i = 2; i <= limit; ++i) {
    if (composite[i]) continue;
    if (i >= kFirstSievingPrime) out.push_back(i);
    for (u64 j = static_cast<u64>(i) * i; j <= limit; j += i) composite[j] = 1;
  }
  return out;
}

u64 isqrt(u64 n) {
  u64 r = static_cast<u64>(std::sqrt(static_cast<double>(n)));
  while (r * r > n) --r;
  while ((r + 1) * (r + 1) <= n) ++r;
  return r;
}

// Upper bound on the n-th prime (Rosser & Schoenfeld / Dusart), n >= 6.
u64 nth_prime_upper_bound(u64 n) {
  if (n < 6) return 13;
  double ln = std::log(static_cast<double>(n));
  double lnln = std::log(ln);
  double bound = n * (ln + lnln - (n >= 688383 ? 0.9484 : 0.0));
  return static_cast<u64>(bound) + 64;
}

class Sieve {
 public:
  // Sieves numbers in [0, limit]. Byte b represents [30b, 30b + 29].
  Sieve(u64 limit, u32 segment_bytes)
      : limit_(limit),
        total_bytes_(limit / 30 + 1),
        segment_bytes_(segment_bytes),
        primes_(sieving_primes(static_cast<u32>(isqrt(limit)))) {}

  u64 total_bytes() const { return total_bytes_; }
  u64 segment_count() const {
    return (total_bytes_ + segment_bytes_ - 1) / segment_bytes_;
  }
  u32 segment_bytes() const { return segment_bytes_; }

  // Sieves segment `index` into `buf` (capacity >= segment_bytes + 8, padded
  // with zeros past the returned length). Returns the segment length in bytes.
  u32 sieve_segment(u64 index, u8* buf) const {
    const Tables& t = tables();
    const u64 lo = index * segment_bytes_;
    const u32 len = static_cast<u32>(
        std::min<u64>(segment_bytes_, total_bytes_ - lo));

    // 1. Pre-sieve by copying the 7*11*13*17 pattern.
    u32 off = static_cast<u32>(lo % kPatternBytes);
    for (u32 done = 0; done < len;) {
      u32 n = std::min(len - done, kPatternBytes - off);
      std::memcpy(buf + done, t.pattern.data() + off, n);
      done += n;
      off = 0;
    }
    std::memset(buf + len, 0, 8);

    // 2. Cross off multiples of each sieving prime p >= 19.
    const u64 lo_num = lo * 30;
    const u64 hi_num = lo_num + static_cast<u64>(len) * 30;  // exclusive
    for (u32 p : primes_) {
      const u64 pp = static_cast<u64>(p) * p;
      if (pp >= hi_num) break;
      // Smallest q >= max(p, ceil(lo_num / p)) with q coprime to 30.
      u64 q = std::max<u64>(p, (lo_num + p - 1) / p);
      int j;
      while ((j = t.residue_index[q % 30]) < 0) ++q;
      cross_off(buf, len, p, q, j, lo);
    }

    // 3. Fix up the edges: 1 is not prime, 7..17 are; clip past limit.
    if (lo == 0) {
      buf[0] &= static_cast<u8>(~1u);  // 1
      for (u32 p : kPreSievePrimes)
        if (p <= limit_) buf[p / 30] |= static_cast<u8>(1u << t.residue_index[p % 30]);
    }
    if (lo + len == total_bytes_) {
      const u64 base = (total_bytes_ - 1) * 30;
      for (int k = 0; k < 8; ++k)
        if (base + kResidues[k] > limit_)
          buf[len - 1] &= static_cast<u8>(~(1u << k));
    }
    return len;
  }

 private:
  static inline void cross_off(u8* buf, u32 len, u32 p, u64 q, int j, u64 lo) {
    const Tables& t = tables();
    const int i = t.residue_index[p % 30];
    const u64 a = p / 30;
    const u64 start = (static_cast<u64>(p) * q) / 30;
    if (start >= lo + len) return;
    u64 pos = start - lo;

    // Advance to the start of a wheel turn (q = 1 mod 30).
    while (j != 0) {
      if (pos >= len) return;
      buf[pos] &= t.mask[i][j];
      pos += a * kGaps[j] + t.carry[i][j];
      j = (j + 1) & 7;
    }

    // Byte offsets of the 8 strikes within one wheel turn; a turn spans p bytes.
    u64 o[8];
    o[0] = 0;
    for (int k = 0; k < 7; ++k) o[k + 1] = o[k] + a * kGaps[k] + t.carry[i][k];
    const u8 m0 = t.mask[i][0], m1 = t.mask[i][1], m2 = t.mask[i][2],
             m3 = t.mask[i][3], m4 = t.mask[i][4], m5 = t.mask[i][5],
             m6 = t.mask[i][6], m7 = t.mask[i][7];
    const u64 o1 = o[1], o2 = o[2], o3 = o[3], o4 = o[4], o5 = o[5],
              o6 = o[6], o7 = o[7];

    if (len > o7) {
      u8* s = buf + pos;
      u8* const end = buf + (len - o7);
      for (; s < end; s += p) {
        s[0] &= m0;
        s[o1] &= m1;
        s[o2] &= m2;
        s[o3] &= m3;
        s[o4] &= m4;
        s[o5] &= m5;
        s[o6] &= m6;
        s[o7] &= m7;
      }
      pos = static_cast<u64>(s - buf);
    }

    // Tail of the last partial turn.
    for (j = 0; pos < len; j = (j + 1) & 7) {
      buf[pos] &= t.mask[i][j];
      pos += a * kGaps[j] + t.carry[i][j];
    }
  }

  u64 limit_;
  u64 total_bytes_;
  u32 segment_bytes_;
  std::vector<u32> primes_;
};

u64 popcount_bytes(const u8* buf, u32 len) {
  u64 count = 0;
  u32 words = (len + 7) / 8;  // buffer is zero-padded past len
  for (u32 w = 0; w < words; ++w) {
    u64 x;
    std::memcpy(&x, buf + 8 * w, 8);
    count += static_cast<u64>(__builtin_popcountll(x));
  }
  return count;
}

// Calls f(prime) for every set bit in buf, in increasing order.
template <typename F>
void for_each_prime(const u8* buf, u32 len, u64 lo_byte, F&& f) {
  for (u32 b = 0; b < len; ++b) {
    u32 bits = buf[b];
    const u64 base = (lo_byte + b) * 30;
    while (bits) {
      int k = __builtin_ctz(bits);
      bits &= bits - 1;
      f(base + kResidues[k]);
    }
  }
}

// Appends n and '\n' to out.
inline char* write_u64(char* out, u64 n) {
  char tmp[24];
  int len = 0;
  do {
    tmp[len++] = static_cast<char>('0' + n % 10);
    n /= 10;
  } while (n);
  while (len) *out++ = tmp[--len];
  *out++ = '\n';
  return out;
}

template <typename Work>
void run_parallel(unsigned threads, Work&& work) {
  std::vector<std::thread> pool;
  for (unsigned t = 1; t < threads; ++t) pool.emplace_back(work);
  work();
  for (auto& th : pool) th.join();
}

// Counts primes in every segment, returns per-segment counts.
std::vector<u32> count_segments(const Sieve& sieve, unsigned threads) {
  const u64 nseg = sieve.segment_count();
  std::vector<u32> counts(nseg);
  std::atomic<u64> next{0};
  run_parallel(threads, [&] {
    std::vector<u8> buf(sieve.segment_bytes() + 8);
    for (u64 s; (s = next.fetch_add(1, std::memory_order_relaxed)) < nseg;) {
      u32 len = sieve.sieve_segment(s, buf.data());
      counts[s] = static_cast<u32>(popcount_bytes(buf.data(), len));
    }
  });
  return counts;
}

// Writes every prime <= limit to `out`, segments sieved in parallel and
// written in order.
void print_primes(u64 limit, u32 segment_bytes, unsigned threads, FILE* out) {
  static const char kSmall[] = "2\n3\n5\n";
  if (limit >= 5) std::fwrite(kSmall, 1, 6, out);
  else if (limit >= 3) std::fwrite(kSmall, 1, 4, out);
  else if (limit >= 2) std::fwrite(kSmall, 1, 2, out);
  if (limit < 7) return;

  Sieve sieve(limit, segment_bytes);
  const u64 nseg = sieve.segment_count();
  // At most 8 primes per byte, at most 21 chars per prime.
  const size_t max_text = static_cast<size_t>(segment_bytes) * 8 * 21;
  std::vector<std::vector<u8>> bufs(threads, std::vector<u8>(segment_bytes + 8));
  std::vector<std::vector<char>> text(threads, std::vector<char>(max_text));
  std::vector<size_t> text_len(threads);

  for (u64 batch = 0; batch < nseg; batch += threads) {
    const unsigned n = static_cast<unsigned>(std::min<u64>(threads, nseg - batch));
    std::vector<std::thread> pool;
    auto job = [&](unsigned t) {
      u64 s = batch + t;
      u32 len = sieve.sieve_segment(s, bufs[t].data());
      char* p = text[t].data();
      for_each_prime(bufs[t].data(), len, s * sieve.segment_bytes(),
                     [&](u64 prime) { p = write_u64(p, prime); });
      text_len[t] = static_cast<size_t>(p - text[t].data());
    };
    for (unsigned t = 1; t < n; ++t) pool.emplace_back(job, t);
    job(0);
    for (auto& th : pool) th.join();
    for (unsigned t = 0; t < n; ++t) std::fwrite(text[t].data(), 1, text_len[t], out);
  }
}

}  // namespace

int main(int argc, char** argv) {
  u64 n = 1000000000ull;
  unsigned threads = std::max(1u, std::thread::hardware_concurrency());
  u32 segment_kb = 256;
  const char* print_path = nullptr;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    auto need = [&](const char* name) -> const char* {
      if (a + 1 >= argc) {
        std::fprintf(stderr, "missing value for %s\n", name);
        std::exit(1);
      }
      return argv[++a];
    };
    if (arg == "--threads") threads = std::max(1, std::atoi(need("--threads")));
    else if (arg == "--segment") segment_kb = std::max(1, std::atoi(need("--segment")));
    else if (arg == "--print") print_path = need("--print");
    else if (arg == "-h" || arg == "--help") {
      std::printf("usage: %s [N] [--threads T] [--segment KB] [--print FILE]\n", argv[0]);
      return 0;
    } else {
      n = std::strtoull(arg.c_str(), nullptr, 10);
      if (n == 0) {
        std::fprintf(stderr, "N must be a positive integer\n");
        return 1;
      }
    }
  }
  const u32 segment_bytes = segment_kb * 1024u;

  auto t0 = std::chrono::steady_clock::now();

  u64 nth = 0;
  constexpr u64 kSmall[3] = {2, 3, 5};
  if (n <= 3) {
    nth = kSmall[n - 1];
  } else {
    const u64 limit = nth_prime_upper_bound(n);
    Sieve sieve(limit, segment_bytes);
    std::vector<u32> counts = count_segments(sieve, threads);

    // Locate the segment containing the n-th prime, then scan it.
    u64 remaining = n - 3;  // 2, 3, 5 are not in the wheel
    u64 seg = 0;
    while (counts[seg] < remaining) remaining -= counts[seg++];
    std::vector<u8> buf(segment_bytes + 8);
    u32 len = sieve.sieve_segment(seg, buf.data());
    for_each_prime(buf.data(), len, seg * segment_bytes, [&](u64 p) {
      if (--remaining == 0) nth = p;
    });
  }

  auto t1 = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(t1 - t0).count();
  std::fprintf(stderr, "prime #%llu = %llu  (%.3f s, %u threads)\n",
               static_cast<unsigned long long>(n),
               static_cast<unsigned long long>(nth), secs, threads);

  if (print_path) {
    FILE* out = std::strcmp(print_path, "-") == 0 ? stdout : std::fopen(print_path, "wb");
    if (!out) {
      std::perror(print_path);
      return 1;
    }
    static char iobuf[1 << 22];
    std::setvbuf(out, iobuf, _IOFBF, sizeof iobuf);
    print_primes(nth, segment_bytes, threads, out);
    if (out != stdout) std::fclose(out);
    auto t2 = std::chrono::steady_clock::now();
    std::fprintf(stderr, "wrote %llu primes to %s  (%.3f s)\n",
                 static_cast<unsigned long long>(n), print_path,
                 std::chrono::duration<double>(t2 - t1).count());
  }
  return 0;
}
