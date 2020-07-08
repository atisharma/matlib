[![PRsWelcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![](https://cljdoc.org/badge/matlib)](https://cljdoc.org/jump/release/matlib)

# matlib

Matlib is a Clojure library of optimisation and control theory tools and
convenience functions based on Neanderthal.

## Motivation

Clojure's REPL-driven workflow is well suited to exploring and manipulating data.
Like MATLAB, this approach should work well for matrix-based numerical computation.

Neanderthal fills a need for a performant matrix library in Clojure, but it is
essentially a thin wrapper around LAPACK functions. This library aims to
furnish Clojure with some higher-level functions and applications, including
system identification, control theory and optimisation tools, without replicating
things already available in Neanderthal.

Why not core.matrix? It is of course subjective. The philosophy of core.matrix
and Neanderthal is different and they serve different needs. I preferred being
close to LAPACK -- the cockroach of numerical computing -- and Neanderthal's
syntax made more sense to me.


## Installing

TBD


## Contributing

Pull requests and bug reports are welcome. 

The code is written in a style that stays close to the mathematics in the
referenced papers where possible. This leads to extensive use of `let`.


### Finished

- Various linear algebra functions like pseudo-inverse, kernel, subspace projections etc.
- (optimisation) L-BFGS and gradient descent
- (system identification) N4SID first, second (biased), and robust algorithms (untested)
- Basis state-space representation, discrete-time integration
- Gramians, Lyapunov equations
- Some convenience functions


### Unfinished

- There are currently **no tests**
- Continuous-time integration of state-space models
- Kroneker product
- MOESP B, D, covariances
- Other system ID algorithms
- Continuous-time state space integration
- Riccati equation solver (but see [this issue](https://github.com/uncomplicate/neanderthal/issues/93))
- Complex matrices (lacking in Neanderthal)

Matlib uses only Neanderthal's native type. It shouldn't be *that* hard to use
Neanderthal's GPU types, but I have not attempted it.


## License

Copyright © 2020 A S Sharma

This program and the accompanying materials are made available under the
terms of the Eclipse Public License 2.0 which is available at
http://www.eclipse.org/legal/epl-2.0.

This Source Code may also be made available under the following Secondary
Licenses when the conditions for such availability set forth in the Eclipse
Public License, v. 2.0 are satisfied: GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at your
option) any later version, with the GNU Classpath Exception which is available
at https://www.gnu.org/software/classpath/license.html.
