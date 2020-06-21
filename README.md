# matlib

Matlib is a Clojure library of convenience matrix functions based on Neanderthal.

I could have called it Mousterian, but I prefer descriptive names.

## Motivation

Clojure's REPL-driven workflow is well suited to exploring and manipulating data.
Like MATLAB, this approach should work well for matrix-based numerical computation.

Neanderthal fills a need for a performant matrix library in Clojure, but it is
essentially a thin wrapper around LAPACK functions. This library aims to
furnish Clojure with some higher-level functions and applications, including
system identification, control theory and optimisation tools.

Why not core.matrix? It is of course subjective. The philosophy of core.matrix
and Neanderthal is different and they serve different needs. I preferred being
close to LAPACK -- the cockroach of numerical computing -- and Neanderthal's
syntax made more sense to me. If it's possible to use well-tested FORTRAN
libraries for numerical computing, I will.

## Contributing

Pull requests and bug reports are welcome. 

The code is written in a style that is meant to read like the mathematics where
possible. This leads to extensive use of `let`, which may not be pretty
Clojure, but is easier to relate back to papers describing the mathematics.


### Finished

- SVD-based linear algebra functions like pseudo-inverse, subspace projections etc.
- N4SID second algorithm
- Basis state-space representation, discrete-time integration
- Some convenience functions


### Unfinished

- There are currently no tests
- Schur decomposition, Kroneker product
- Other system ID algorithms
- Continuous-time state space integration
- Control theory tools such as a Riccati equation solver, Lyapunov solver
- optimisers such as an L-BFGS implementation
- Saving and loading of matrices
- Complex matrices (lacking in Neanderthal)

Matlib uses only Neanderthal's native type. It shouldn't be *that* hard to use
Neanderthal's GPU types, but I have not attempted the refactoring and testing
necessary to successfully support it.


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
