#pragma once
#include <array>
#include <cstddef>
#include <type_traits>

namespace tinyfit::pso {

/* Header-only Particle Swarm Optimization.
 *   T           : scalar type (float/double)
 *   Dim         : search space dimension
 *   Particles   : swarm size
 *   Iterations  : iteration count
 *   Objective   : callable  T(const std::array<T,Dim>&)
 *   RNG         : object with double operator()() in [0,1)
 */
template <typename T, std::size_t Dim, std::size_t Particles,
          std::size_t Iterations, auto Objective, typename RNG>
[[nodiscard]] constexpr auto run(RNG rng) {
  static_assert(std::is_arithmetic_v<T>);
  constexpr T w = T(0.729844);  // Clerc–Kennedy χ
  constexpr T c1 = T(1.49618);
  constexpr T c2 = T(1.49618);

  using Vec = std::array<T, Dim>;
  std::array<Vec, Particles> pos{}, vel{}, pbest{};
  std::array<T, Particles> pbest_val{};
  Vec gbest{};
  T gbest_val{};

  auto rand_range = [&](double r) {
    return T(r * (rng() * 2.0 - 1.0));
  };

  /* --- initialise swarm --- */
  for (std::size_t p = 0; p < Particles; ++p) {
    for (std::size_t d = 0; d < Dim; ++d) {
      pos[p][d] = rand_range(5.0);
      vel[p][d] = rand_range(1.0);
    }
    pbest[p] = pos[p];
    pbest_val[p] = Objective(pos[p]);
    if (p == 0 || pbest_val[p] < gbest_val) {
      gbest_val = pbest_val[p];
      gbest = pbest[p];
    }
  }

  /* --- main loop --- */
  for (std::size_t it = 0; it < Iterations; ++it) {
    for (std::size_t p = 0; p < Particles; ++p) {
      for (std::size_t d = 0; d < Dim; ++d) {
        double r1 = rng(), r2 = rng();
        vel[p][d] = w * vel[p][d] + c1 * r1 * (pbest[p][d] - pos[p][d]) +
                    c2 * r2 * (gbest[d] - pos[p][d]);
        pos[p][d] += vel[p][d];
      }
      T val = Objective(pos[p]);
      if (val < pbest_val[p]) {
        pbest_val[p] = val;
        pbest[p] = pos[p];
        if (val < gbest_val) {
          gbest_val = val;
          gbest = pos[p];
        }
      }
    }
  }
  return std::pair{gbest, gbest_val};
}

}  // namespace tinyfit::pso
