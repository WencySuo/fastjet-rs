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

use criterion::{Criterion, criterion_group, criterion_main};

fn criterion_benchmark(c: &mut Criterion) {
    //TODO: prolly optimize reading lol
    //
    c.bench_function("example_01", |b| {
        b.iter(|| {
            let file = File::open("./examples/data/single-event.dat").unwrap();
            let mut reader = BufReader::new(file);

            let mut input_particles = Vec::new();
            let mut line = String::new();

            while reader.read_line(&mut line).unwrap() != 0 {
                let mut it = line.split_whitespace();

                let px: f64 = it.next().unwrap().parse().unwrap();
                let py: f64 = it.next().unwrap().parse().unwrap();
                let pz: f64 = it.next().unwrap().parse().unwrap();
                let e: f64 = it.next().unwrap().parse().unwrap();

                input_particles.push(PseudoJet::new(px, py, pz, e));

                line.clear();
            }

            // let mut times: Vec<u128> = Vec::new();

            let stdout = stdout();
            let mut out = BufWriter::new(stdout.lock());

            // let now = std::time::Instant::now();

            //TODO: more default constructors

            //create jet def
            let _r: f64 = 0.6;
            let jet_def = JetDefinition::new(
                Algorithm::AntiKt,
                0.6,
                RecombinationScheme::EScheme,
                Strategy::N2Plain,
                None,
            );

            let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

            clust_seq.initialize_and_run_no_decant();

            let pmin = 5.0;

            let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

            let inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

            //let elapsed = now.elapsed().as_nanos();
            //times.push(elapsed);

            writeln!(
                out,
                "{:>5} {:>15} {:>15} {:>15}",
                "jet #", "rapiddity", "phi", "pt"
            )
            .unwrap();

            for (i, jet) in inclusive_jets.iter().enumerate() {
                writeln!(
                    out,
                    "{:>5} {:>15.8} {:>15.8} {:>15.8}",
                    i,
                    jet.rap(),
                    jet.phi(),
                    jet.pt()
                )
                .unwrap();
            }
            // let elapsed = now.elapsed()?.as_nanos();
            // times.push(elapsed);
            //}

            // let mean = times.iter().sum::<u128>() as f64 / times.len() as f64;

            // let variance_sum: f64 = times
            //     .iter()
            //     .map(|&x| {
            //         let diff = x as f64 - mean;
            //         diff * diff
            //     })
            //     .sum();

            // let sample_variance = variance_sum as f64 / (times.len() - 1) as f64;

            // let uncertainty = sample_variance.sqrt() / (times.len() as f64).sqrt();

            // writeln!(
            //     out,
            //     "The average time for 1000 iterations is {:.3} µs ± {:.3} µs",
            //     mean / 1000.0,
            //     uncertainty / 1000.0
            // )?;
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
