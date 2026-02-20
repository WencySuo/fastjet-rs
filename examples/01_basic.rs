use fastjet_rs::cluster_sequence::Algorithm;
use fastjet_rs::cluster_sequence::ClusterSequence;
use fastjet_rs::cluster_sequence::JetDefinition;
use fastjet_rs::cluster_sequence::RecombinationScheme;
use fastjet_rs::cluster_sequence::Strategy;
use fastjet_rs::pseudo_jet::PseudoJet;
use std::fs::File;
use std::io::stdout;
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};

use std::time::SystemTime;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //TODO: prolly optimize reading lol
    //
    let now = SystemTime::now();

    let file = File::open("./examples/data/single-event.dat")?;
    let mut reader = BufReader::new(file);

    let mut input_particles = Vec::new();
    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        let mut it = line.split_whitespace();

        let px: f64 = it.next().unwrap().parse().unwrap();
        let py: f64 = it.next().unwrap().parse().unwrap();
        let pz: f64 = it.next().unwrap().parse().unwrap();
        let e: f64 = it.next().unwrap().parse().unwrap();

        input_particles.push(PseudoJet::new(px, py, pz, e));

        line.clear();
    }

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

    let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

    let inclusive_jets: &mut Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

    println!("Time elapse is {:?}", now.elapsed());

    let stdout = stdout();
    let mut out = BufWriter::new(stdout.lock());

    writeln!(
        out,
        "{:>5} {:>15} {:>15} {:>15}",
        "jet #", "rapiddity", "phi", "pt"
    )?;

    for (i, jet) in inclusive_jets.iter_mut().enumerate() {
        writeln!(
            out,
            "{:>5} {:>15.8} {:>15.8} {:>15.8}",
            i,
            jet.rap(),
            jet.phi(),
            jet.pt()
        )?;
    }

    Ok(())
}
