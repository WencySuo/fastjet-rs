### FastJet-rs TODO List

## Input features:
- Implement support for reading input from HepMC files
- Add support for batching
- Add python bindings

## Cluster Sequencing:
- Implement support for different algorithms other than AntiKt
- Implement support for different clustering algorithms (Notable N2Tiled
- Investigate using Rust nightly for 'std::simd' features on floats
- Implement exclusive jets and other functions acting on already clustered jets
- Start playing around with EE jets

## Visualization
- Investigate Rust GUI crates for visualization of jet clustering and output 

## Testing
- Include more example files with different data
- Run CI/CD script to verify correctness with FastJet through integration tests
- backfill unit tests for ClusterSequence, ProxyJet

... and much more (hopefully).
