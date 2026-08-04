// Independent exhaustive verifier for the GF(2) <15,15,17> rank-2256 certificate.
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {
constexpr std::size_t N = 15;
constexpr std::size_t M = 15;
constexpr std::size_t P = 17;
constexpr std::size_t ExpectedRank = 2256;
constexpr std::size_t AB = N * M;
constexpr std::size_t BC = M * P;
constexpr std::size_t AC = N * P;

struct Term {
    cpp_int u;
    cpp_int v;
    cpp_int w;
};

std::string trim(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

cpp_int parse_decimal(const std::string& text) {
    if (text.empty() || !std::all_of(text.begin(), text.end(), [](unsigned char c) { return std::isdigit(c); })) {
        throw std::runtime_error("invalid nonnegative decimal integer: " + text);
    }
    cpp_int value = 0;
    for (const char c : text) {
        value *= 10;
        value += static_cast<unsigned>(c - '0');
    }
    return value;
}

std::vector<std::size_t> bit_positions(const cpp_int& value, std::size_t bits) {
    std::vector<std::size_t> positions;
    for (std::size_t bit = 0; bit < bits; ++bit) {
        if (boost::multiprecision::bit_test(value, bit)) positions.push_back(bit);
    }
    return positions;
}

std::size_t bit_count(cpp_int value) {
    std::size_t count = 0;
    while (value != 0) {
        if ((value & 1) != 0) ++count;
        value >>= 1;
    }
    return count;
}

struct Parsed {
    std::size_t declared_rank = 0;
    std::vector<Term> terms;
    std::size_t duplicate_triples = 0;
};

Parsed parse_certificate(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open certificate: " + path);

    Parsed parsed;
    bool have_rank = false;
    std::unordered_set<std::string> unique;
    std::string raw;
    std::size_t line_number = 0;

    while (std::getline(input, raw)) {
        ++line_number;
        const std::string line = trim(raw);
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::vector<std::string> tokens;
        for (std::string token; fields >> token;) tokens.push_back(token);
        if (!tokens.empty() && tokens.front() == "R") tokens.erase(tokens.begin());

        if (!have_rank) {
            if (tokens.size() != 1) throw std::runtime_error("line " + std::to_string(line_number) + ": expected rank header");
            parsed.declared_rank = static_cast<std::size_t>(std::stoull(tokens[0]));
            have_rank = true;
            continue;
        }
        if (tokens.size() != 3) throw std::runtime_error("line " + std::to_string(line_number) + ": expected U V W");
        const std::string key = tokens[0] + "\n" + tokens[1] + "\n" + tokens[2];
        if (!unique.insert(key).second) ++parsed.duplicate_triples;
        parsed.terms.push_back({parse_decimal(tokens[0]), parse_decimal(tokens[1]), parse_decimal(tokens[2])});
    }
    if (!have_rank) throw std::runtime_error("missing rank header");
    return parsed;
}
}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "matmul_15x15x17_rank2256_gf2.txt";
    const auto started = std::chrono::steady_clock::now();
    try {
        const Parsed parsed = parse_certificate(path);
        const cpp_int limit_u = cpp_int(1) << AB;
        const cpp_int limit_v = cpp_int(1) << BC;
        const cpp_int limit_w = cpp_int(1) << AC;

        std::size_t bad_terms = 0;
        std::vector<cpp_int> residual(AB * BC);
        std::uint64_t expanded_uv_pairs = 0;

        for (const Term& term : parsed.terms) {
            if (term.u <= 0 || term.u >= limit_u || term.v <= 0 || term.v >= limit_v || term.w <= 0 || term.w >= limit_w) {
                ++bad_terms;
                continue;
            }
            const auto u_bits = bit_positions(term.u, AB);
            const auto v_bits = bit_positions(term.v, BC);
            expanded_uv_pairs += static_cast<std::uint64_t>(u_bits.size()) * static_cast<std::uint64_t>(v_bits.size());
            for (const std::size_t a : u_bits) {
                const std::size_t base = a * BC;
                for (const std::size_t b : v_bits) residual[base + b] ^= term.w;
            }
        }

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < M; ++j) {
                const std::size_t a = i * M + j;
                const std::size_t base = a * BC;
                for (std::size_t k = 0; k < P; ++k) {
                    const std::size_t b = j * P + k;
                    residual[base + b] ^= cpp_int(1) << (i * P + k);
                }
            }
        }

        std::size_t residual_pair_slices = 0;
        std::size_t residual_coefficients = 0;
        for (const cpp_int& mask : residual) {
            if (mask != 0) {
                ++residual_pair_slices;
                residual_coefficients += bit_count(mask);
            }
        }

        const bool verified = parsed.declared_rank == ExpectedRank && parsed.terms.size() == ExpectedRank &&
                              bad_terms == 0 && parsed.duplicate_triples == 0 && residual_coefficients == 0;
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();

        std::cout << "{\n"
                  << "  \"certificate\": \"" << path << "\",\n"
                  << "  \"tensor\": [15, 15, 17],\n"
                  << "  \"expected_rank\": " << ExpectedRank << ",\n"
                  << "  \"declared_rank\": " << parsed.declared_rank << ",\n"
                  << "  \"parsed_terms\": " << parsed.terms.size() << ",\n"
                  << "  \"nonzero_and_in_range\": " << (bad_terms == 0 ? "true" : "false") << ",\n"
                  << "  \"bad_terms\": " << bad_terms << ",\n"
                  << "  \"duplicate_triples\": " << parsed.duplicate_triples << ",\n"
                  << "  \"expanded_uv_pairs\": " << expanded_uv_pairs << ",\n"
                  << "  \"target_coefficients\": " << N * M * P << ",\n"
                  << "  \"brent_coefficients\": " << AB * BC * AC << ",\n"
                  << "  \"residual_pair_slices\": " << residual_pair_slices << ",\n"
                  << "  \"residual_coefficients\": " << residual_coefficients << ",\n"
                  << "  \"verified\": " << (verified ? "true" : "false") << ",\n"
                  << "  \"seconds\": " << seconds << "\n"
                  << "}\n";
        return verified ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "verification error: " << error.what() << '\n';
        return 2;
    }
}
