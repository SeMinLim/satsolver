# SAT Solver

- Conflict-Driven Clause Learning
  - First UIP (Unit Implication Point) Clause Learning [Based on Chaff]
- Boolean Constraint Propagation [Based on MiniSAT]
  - Phase Saving (A component caching technique)
- Lazy Data Structure
  - The Watched Literal [Based on Chaff]
- Literal Block Distance Scoring Scheme [Based on Glucose]
- Dynamic Search Restart [Based on Glucose]
- Conflict-Driven Branching Heuristics
  - VSIDS (Variable State Independent Decaying Sum) with Heap Data Structure (Max-Heap) [Based on Chaff]
- Clause Deletion Strategies
  - Reducing based on LBD
- Rephasing [Based on CaDiCaL]

------
### Mandatory Papers

- **[CDCL]** J. Marques-Silva, L. Inês, and M. Sharad, "Conflict-driven clause learning SAT solvers," Handbook of satisfiability, IOS press, 2021, 133-182.
- **[Chaff]** MW. Moskewicz et al, "Chaff: Engineering an efficient SAT solver," Proceedings of the 38th annual Design Automation Conference, 2001.
- **[MiniSAT]** N. Eén and S. Niklas, "An extensible SAT-solver," Lecture notes in computer science 2919.2004 (2004): 502-518.
- **[Phase Saving technique]** K. Pipatsrisawat and D. Adnan, "A lightweight component caching scheme for satisfiability solvers," Theory and Applications of Satisfiability Testing–SAT 2007: 10th International Conference, Lisbon, Portugal, May 28-31, 2007.
- **[Glucose]** G. Audemard and S. Laurent, "Predicting learnt clauses quality in modern SAT solvers," Twenty-first international joint conference on artificial intelligence, 2009.
- **[CaDiCaL]** A. Biere, "Cadical, lingeling, plingeling, treengeling and yalsat entering the sat competition 2017," Proceedings of SAT Competition 14 (2017): 316-336.
