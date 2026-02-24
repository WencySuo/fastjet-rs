# fastjet-rs

Rust implementation of [FastJet](https://fastjet.fr/) originally implemented in C++ for jet finding. This is a work in progress and we are actively working toward feature compatibility with FastJet in a 1.0 release.

Here are the currently supported algorithms for pp collisions: 
- Inclusive and exclusive variants of the Kt jet algorithm 
- Inclusive and exclusive variants of the AntiKt jet algorithm 
- Inclusive and exclusive variants of the Cambridge/Aachen jet algorithm 
- Generalized kt jet algorithm  

## Setup 
Download the library from crates.io:
   ```bash
   cargo add fastjet-rs
   ```

## Basic Usage 
```rust 
use fastjet_rs::*; 

// define vector of input particles 
let input_particles = vec![
    PseudoJet::new(0.1, 3.0, 0.0, 52.0), 
    PseudoJet::new(1.0, 2.0, -4.0, 6.0), 
    PseudoJet::new(-0.1, 82.0, -3.2, 70.0), 
];

// define the jet finding parameters
let jet_def = fastjJetDefinition::new(
    Algorithm::AntiKt,
    0.6,
    RecombinationScheme::EScheme,
    Strategy::N2Plain,
    None, // no need for extra parameters with AntiKt
);

// run the clustering
let mut clust_seq = ClusterSequence::new(input_particles, jet_def);
clust_seq.initialize_and_run_no_decant();

//  define cutoff and obtain output inclusive list of jets 
let pmin = 5.0; 
let mut inclusive_jets = clust_seq.inclusive_jets(pmin);
let inclusive_jets = PseudoJet::sorted_by_pt(&mut inclusive_jets);
```
   
## Roadmap 
### Input features:
- Implement support for reading input from HepMC files
- Add support for batching
- Add python bindings

### Cluster Sequencing:
- Implement support for different algorithms other than AntiKt
- Implement support for different clustering algorithms (Notable N2Tiled
- Investigate using Rust nightly for 'std::simd' features on floats
- Implement exclusive jets and other functions acting on already clustered jets
- Start playing around with EE jets

### Visualization
- Investigate Rust GUI crates for visualization of jet clustering and output 

### Testing
- Include more example files with different data
- Run CI/CD script to verify correctness with FastJet through integration tests
- backfill unit tests for ClusterSequence, ProxyJet

... and much more!
