/**
 * Nelder-Mead Simplex Optimization
 *
 * Based on fmin library by Ben Frederickson (BSD-3-Clause)
 * https://github.com/benfred/fmin
 *
 * Patched to fix variable scope issue with 'i' at line 132 of original.
 */

/**
 * Compute dot product of two vectors.
 */
function dot(a, b) {
  let ret = 0;
  for (let idx = 0; idx < a.length; ++idx) {
    ret += a[idx] * b[idx];
  }
  return ret;
}

/**
 * Compute weighted sum of two vectors: out = c1*a + c2*b
 */
function weightedSum(out, c1, a, c2, b) {
  for (let idx = 0; idx < out.length; ++idx) {
    out[idx] = c1 * a[idx] + c2 * b[idx];
  }
}

/**
 * Minimizes a function using the Nelder-Mead downhill simplex method.
 *
 * @param {Function} f - Function to minimize, takes array and returns scalar
 * @param {Array<number>} x0 - Initial parameter values
 * @param {Object} [parameters] - Optional parameters
 * @param {number} [parameters.maxIterations] - Maximum iterations (default: 200 * dimensions)
 * @param {number} [parameters.minErrorDelta] - Convergence tolerance on function value
 * @param {number} [parameters.minTolerance] - Convergence tolerance on parameter values
 * @returns {Object} { fx: minimum value, x: minimizing parameters }
 */
export function nelderMead(f, x0, parameters) {
  parameters = parameters || {};

  const maxIterations = parameters.maxIterations || x0.length * 200;
  const nonZeroDelta = parameters.nonZeroDelta || 1.05;
  const zeroDelta = parameters.zeroDelta || 0.001;
  const minErrorDelta = parameters.minErrorDelta || 1e-6;
  const minTolerance = parameters.minTolerance || 1e-5;
  const rho = parameters.rho !== undefined ? parameters.rho : 1;
  const chi = parameters.chi !== undefined ? parameters.chi : 2;
  const psi = parameters.psi !== undefined ? parameters.psi : -0.5;
  const sigma = parameters.sigma !== undefined ? parameters.sigma : 0.5;
  let maxDiff;

  // Initialize simplex
  const N = x0.length;
  const simplex = new Array(N + 1);
  simplex[0] = x0.slice();
  simplex[0].fx = f(x0);
  simplex[0].id = 0;

  for (let idx = 0; idx < N; ++idx) {
    const point = x0.slice();
    point[idx] = point[idx] ? point[idx] * nonZeroDelta : zeroDelta;
    simplex[idx + 1] = point;
    simplex[idx + 1].fx = f(point);
    simplex[idx + 1].id = idx + 1;
  }

  function updateSimplex(value) {
    for (let idx = 0; idx < value.length; idx++) {
      simplex[N][idx] = value[idx];
    }
    simplex[N].fx = value.fx;
  }

  const sortOrder = (a, b) => a.fx - b.fx;

  const centroid = x0.slice();
  const reflected = x0.slice();
  const contracted = x0.slice();
  const expanded = x0.slice();

  for (let iteration = 0; iteration < maxIterations; ++iteration) {
    simplex.sort(sortOrder);

    if (parameters.history) {
      const sortedSimplex = simplex.map((x) => {
        const state = x.slice();
        state.fx = x.fx;
        state.id = x.id;
        return state;
      });
      sortedSimplex.sort((a, b) => a.id - b.id);

      parameters.history.push({
        x: simplex[0].slice(),
        fx: simplex[0].fx,
        simplex: sortedSimplex,
      });
    }

    maxDiff = 0;
    for (let idx = 0; idx < N; ++idx) {
      maxDiff = Math.max(maxDiff, Math.abs(simplex[0][idx] - simplex[1][idx]));
    }

    if (Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta && maxDiff < minTolerance) {
      break;
    }

    // Compute the centroid of all but the worst point in the simplex
    for (let idx = 0; idx < N; ++idx) {
      centroid[idx] = 0;
      for (let jdx = 0; jdx < N; ++jdx) {
        centroid[idx] += simplex[jdx][idx];
      }
      centroid[idx] /= N;
    }

    // Reflect the worst point past the centroid
    const worst = simplex[N];
    weightedSum(reflected, 1 + rho, centroid, -rho, worst);
    reflected.fx = f(reflected);

    // If the reflected point is the best seen, possibly expand
    if (reflected.fx < simplex[0].fx) {
      weightedSum(expanded, 1 + chi, centroid, -chi, worst);
      expanded.fx = f(expanded);
      if (expanded.fx < reflected.fx) {
        updateSimplex(expanded);
      } else {
        updateSimplex(reflected);
      }
    }
    // If the reflected point is worse than the second worst, contract
    else if (reflected.fx >= simplex[N - 1].fx) {
      let shouldReduce = false;

      if (reflected.fx > worst.fx) {
        // Inside contraction
        weightedSum(contracted, 1 + psi, centroid, -psi, worst);
        contracted.fx = f(contracted);
        if (contracted.fx < worst.fx) {
          updateSimplex(contracted);
        } else {
          shouldReduce = true;
        }
      } else {
        // Outside contraction
        weightedSum(contracted, 1 - psi * rho, centroid, psi * rho, worst);
        contracted.fx = f(contracted);
        if (contracted.fx < reflected.fx) {
          updateSimplex(contracted);
        } else {
          shouldReduce = true;
        }
      }

      if (shouldReduce) {
        // If we don't contract here, we're done
        if (sigma >= 1) break;

        // Do a reduction - FIX: properly declare loop variable
        for (let sIdx = 1; sIdx < simplex.length; ++sIdx) {
          weightedSum(simplex[sIdx], 1 - sigma, simplex[0], sigma, simplex[sIdx]);
          simplex[sIdx].fx = f(simplex[sIdx]);
        }
      }
    } else {
      updateSimplex(reflected);
    }
  }

  simplex.sort(sortOrder);
  return { fx: simplex[0].fx, x: simplex[0] };
}
