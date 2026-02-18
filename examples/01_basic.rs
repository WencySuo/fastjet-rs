use fastjet_rs::cluster_sequence::Algorithm;
use fastjet_rs::cluster_sequence::ClusterSequence;
use fastjet_rs::cluster_sequence::JetDefinition;
use fastjet_rs::cluster_sequence::RecombinationScheme;
use fastjet_rs::cluster_sequence::Strategy;
use fastjet_rs::pseudo_jet::PseudoJet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //TODO: prolly optimize reading lol

    // read in input particles
    let file_path = Path::new("./examples/data/single-event.dat");
    let file = File::open(&file_path)?;
    let reader = BufReader::new(file);

    let mut input_particles: Vec<PseudoJet> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        //create pseudojet
        let floats = line.trim().split_whitespace();
        let floats = floats.map(|s| s.parse().unwrap()).collect::<Vec<f64>>();

        let px = floats[0];
        let py = floats[1];
        let pz = floats[2];
        let e = floats[3];

        //println!("Particle: ({}, {}, {}, {})", px, py, pz, e);
        input_particles.push(PseudoJet::new(px, py, pz, e));
    }

    println!("Input Particles size {}", input_particles.len());

    //TODO: more default constructors

    //create jet def
    let _r: f64 = 0.6;
    let jet_def = JetDefinition::new(
        Algorithm::AntiKt,
        0.6,
        RecombinationScheme::EScheme,
        Strategy::N2Plain,
    );

    let mut clust_seq = ClusterSequence::new(input_particles, jet_def);

    clust_seq.initialize_and_run_no_decant();

    let pmin = 5.0;

    println!("Particles size {}", clust_seq.particles.len());

    let mut inclusive_jets: Vec<PseudoJet> = clust_seq.inclusive_jets(pmin);

    println!("Inclusive Jets size {}", inclusive_jets.len());
    let inclusive_jets: &mut Vec<PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

    println!(
        "{:>5} {:>15} {:>15} {:>15}\n",
        "jet #", "rapidity", "phi", "pt"
    );

    for (i, jet) in inclusive_jets.iter_mut().enumerate() {
        println!(
            "{:>5} {:>15.8} {:>15.8} {:>15.8}",
            i,
            jet.rap(),
            jet.phi(),
            jet.pt()
        );
    }

    Ok(())
}
