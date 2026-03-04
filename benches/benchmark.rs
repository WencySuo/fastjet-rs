use criterion::BenchmarkId;
use fastjet_rs::cluster_sequence::Algorithm;
use fastjet_rs::cluster_sequence::ClusterSequence;
use fastjet_rs::cluster_sequence::JetDefinition;
use fastjet_rs::cluster_sequence::RecombinationScheme;
use fastjet_rs::cluster_sequence::Strategy;
use fastjet_rs::pseudo_jet::PseudoJet;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor};
use std::time::Duration;
// use std::io::{BufWriter, Write};

use criterion::{Criterion, criterion_group, criterion_main};

fn count_exact_newlines(s: &str) -> usize {
    s.chars().filter(|&c| c == '\n').count()
}

fn run_clustering(event: &String, strategy: Strategy) {
    let event_cursor = Cursor::new(event.as_bytes());
    //let file = File::open("./examples/data/single-event.dat").unwrap();
    let mut reader = BufReader::new(event_cursor);

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

    // let stdout = stdout();
    // let mut out = BufWriter::new(stdout.lock());

    // let now = std::time::Instant::now();

    //TODO: more default constructors

    //create jet def
    let _r: f64 = 0.6;
    let jet_def = JetDefinition::new(
        Algorithm::AntiKt,
        0.6,
        RecombinationScheme::EScheme,
        strategy,
        None,
    );

    let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

    clust_seq.initialize_and_run_no_decant();

    let pmin = 5.0;

    let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

    let _inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);
}

fn example_bench(c: &mut Criterion) {
    c.bench_function("example_01_simd", |b| {
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

            // let stdout = stdout();
            // let mut out = BufWriter::new(stdout.lock());

            // let now = std::time::Instant::now();

            //TODO: more default constructors

            //create jet def
            let _r: f64 = 0.6;
            let jet_def = JetDefinition::new(
                Algorithm::AntiKt,
                0.6,
                RecombinationScheme::EScheme,
                Strategy::N2Simd,
                None,
            );

            let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

            clust_seq.initialize_and_run_no_decant();

            let pmin = 5.0;

            let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

            let _inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

            //let elapsed = now.elapsed().as_nanos();
            //times.push(elapsed);

            // writeln!(
            //     out,
            //     "{:>5} {:>15} {:>15} {:>15}",
            //     "jet #", "rapiddity", "phi", "pt"
            // )
            // .unwrap();

            // for (i, jet) in inclusive_jets.iter().enumerate() {
            //     writeln!(
            //         out,
            //         "{:>5} {:>15.8} {:>15.8} {:>15.8}",
            //         i,
            //         jet.rap(),
            //         jet.phi(),
            //         jet.pt()
            //     )
            //     .unwrap();
            // }
            // let elapsed = now.elapsed()?.as_nanos();
            // times.push(elapsed);
            //}
        })
    });

    c.bench_function("example_01_tiled", |b| {
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

            // let stdout = stdout();
            // let mut out = BufWriter::new(stdout.lock());

            // let now = std::time::Instant::now();

            //TODO: more default constructors

            //create jet def
            let _r: f64 = 0.6;
            let jet_def = JetDefinition::new(
                Algorithm::AntiKt,
                0.6,
                RecombinationScheme::EScheme,
                Strategy::N2Tiling,
                None,
            );

            let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

            clust_seq.initialize_and_run_no_decant();

            let pmin = 5.0;

            let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

            let _inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

            //let elapsed = now.elapsed().as_nanos();
            //times.push(elapsed);

            // writeln!(
            //     out,
            //     "{:>5} {:>15} {:>15} {:>15}",
            //     "jet #", "rapiddity", "phi", "pt"
            // )
            // .unwrap();

            // for (i, jet) in inclusive_jets.iter().enumerate() {
            //     writeln!(
            //         out,
            //         "{:>5} {:>15.8} {:>15.8} {:>15.8}",
            //         i,
            //         jet.rap(),
            //         jet.phi(),
            //         jet.pt()
            //     )
            //     .unwrap();
            // }
            // let elapsed = now.elapsed()?.as_nanos();
            // times.push(elapsed);
            //}
        })
    });

    c.bench_function("example_01_simple", |b| {
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

            // let stdout = stdout();
            // let mut out = BufWriter::new(stdout.lock());

            // let now = std::time::Instant::now();

            //TODO: more default constructors

            //create jet def
            let _r: f64 = 0.6;
            let jet_def = JetDefinition::new(
                Algorithm::AntiKt,
                0.6,
                RecombinationScheme::EScheme,
                Strategy::N2Simd,
                None,
            );

            let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

            clust_seq.initialize_and_run_no_decant();

            let pmin = 5.0;

            let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

            let _inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

            //let elapsed = now.elapsed().as_nanos();
            //times.push(elapsed);

            // writeln!(
            //     out,
            //     "{:>5} {:>15} {:>15} {:>15}",
            //     "jet #", "rapiddity", "phi", "pt"
            // )
            // .unwrap();

            // for (i, jet) in inclusive_jets.iter().enumerate() {
            //     writeln!(
            //         out,
            //         "{:>5} {:>15.8} {:>15.8} {:>15.8}",
            //         i,
            //         jet.rap(),
            //         jet.phi(),
            //         jet.pt()
            //     )
            //     .unwrap();
            // }
            // let elapsed = now.elapsed()?.as_nanos();
            // times.push(elapsed);
            //}
        })
    });

    //         // let mean = times.iter().sum::<u128>() as f64 / times.len() as f64;

    //         // let variance_sum: f64 = times
    //         //     .iter()
    //         //     .map(|&x| {
    //         //         let diff = x as f64 - mean;
    //         //         diff * diff
    //         //     })
    //         //     .sum();

    //         // let sample_variance = variance_sum as f64 / (times.len() - 1) as f64;

    //         // let uncertainty = sample_variance.sqrt() / (times.len() as f64).sqrt();

    //         // writeln!(
    //         //     out,
    //         //     "The average time for 1000 iterations is {:.3} µs ± {:.3} µs",
    //         //     mean / 1000.0,
    //         //     uncertainty / 1000.0
    //         // )?;
    //     })
    // });
}

fn bootstrap_bench(c: &mut Criterion) {
    //instead of loading the file we will pass to bench function
    // each string as a File stream to ensure fair reading of files

    // the file will have 100 events, we guarentee that there is an event
    // in every bin of width 10 from 0 to 1000 thus properly testing all events evenly
    let file = File::open("./examples/data/benchmark.dat").unwrap();
    let reader = BufReader::new(file);
    let mut events: Vec<String> = Vec::new();
    let mut current_event = String::new();

    for line in reader.lines() {
        let line = line.unwrap();
        let trimmed = line.trim();

        if trimmed == "#END" {
            events.push(current_event.clone());
            current_event.clear();
            continue;
        }

        current_event.push_str(trimmed);
        current_event.push('\n');
    }

    let mut group_simple = c.benchmark_group("Simple_N2_cluster");

    // let groups = vec![
    //     (group_simple, Strategy::N2Plain),
    //     (group_simd, Strategy::N2Simd),
    //     (group_tiled, Strategy::N2Tiling),
    // ];

    for event in &events {
        let strategy = Strategy::N2Plain;
        group_simple.throughput(criterion::Throughput::Elements(
            count_exact_newlines(event) as u64
        ));
        group_simple.measurement_time(Duration::from_secs(30));
        group_simple.sample_size(500);
        group_simple.bench_with_input(
            BenchmarkId::from_parameter(count_exact_newlines(event) as u64),
            event,
            |b, event| b.iter(|| run_clustering(event, strategy)),
        );
    }

    group_simple.finish();

    let mut group_simd = c.benchmark_group("N2_Cluster_SIMD");

    for event in &events {
        let strategy = Strategy::N2Simd;
        group_simd.throughput(criterion::Throughput::Elements(
            count_exact_newlines(event) as u64
        ));
        group_simd.measurement_time(Duration::from_secs(30));
        group_simd.sample_size(500);
        group_simd.bench_with_input(
            BenchmarkId::from_parameter(count_exact_newlines(event) as u64),
            event,
            |b, event| b.iter(|| run_clustering(event, strategy)),
        );
    }

    group_simd.finish();

    let mut group_tiled = c.benchmark_group("Tiled_N2_Cluster");

    for event in &events {
        let strategy = Strategy::N2Tiling;
        group_tiled.throughput(criterion::Throughput::Elements(
            count_exact_newlines(event) as u64
        ));
        group_tiled.measurement_time(Duration::from_secs(30));
        group_tiled.sample_size(500);
        group_tiled.bench_with_input(
            BenchmarkId::from_parameter(count_exact_newlines(event) as u64),
            event,
            |b, event| b.iter(|| run_clustering(event, strategy)),
        );
    }

    group_tiled.finish();

    // c.bench_function("example_01", |b| {
    //     b.iter(|| {
    //         let file = File::open("./examples/data/single-event.dat").unwrap();
    //         let mut reader = BufReader::new(file);

    //         let mut input_particles = Vec::new();
    //         let mut line = String::new();

    //         while reader.read_line(&mut line).unwrap() != 0 {
    //             let mut it = line.split_whitespace();

    //             let px: f64 = it.next().unwrap().parse().unwrap();
    //             let py: f64 = it.next().unwrap().parse().unwrap();
    //             let pz: f64 = it.next().unwrap().parse().unwrap();
    //             let e: f64 = it.next().unwrap().parse().unwrap();

    //             input_particles.push(PseudoJet::new(px, py, pz, e));

    //             line.clear();
    //         }

    //         // let mut times: Vec<u128> = Vec::new();

    //         // let stdout = stdout();
    //         // let mut out = BufWriter::new(stdout.lock());

    //         // let now = std::time::Instant::now();

    //         //TODO: more default constructors

    //         //create jet def
    //         let _r: f64 = 0.6;
    //         let jet_def = JetDefinition::new(
    //             Algorithm::AntiKt,
    //             0.6,
    //             RecombinationScheme::EScheme,
    //             Strategy::N2Plain,
    //             None,
    //         );

    //         let mut clust_seq = ClusterSequence::new(input_particles.clone(), jet_def);

    //         clust_seq.initialize_and_run_no_decant();

    //         let pmin = 5.0;

    //         let mut inclusive_jets: Vec<&PseudoJet> = clust_seq.inclusive_jets(pmin);

    //         let _inclusive_jets: &Vec<&PseudoJet> = PseudoJet::sorted_by_pt(&mut inclusive_jets);

    //         //let elapsed = now.elapsed().as_nanos();
    //         //times.push(elapsed);

    //         // writeln!(
    //         //     out,
    //         //     "{:>5} {:>15} {:>15} {:>15}",
    //         //     "jet #", "rapiddity", "phi", "pt"
    //         // )
    //         // .unwrap();

    //         // for (i, jet) in inclusive_jets.iter().enumerate() {
    //         //     writeln!(
    //         //         out,
    //         //         "{:>5} {:>15.8} {:>15.8} {:>15.8}",
    //         //         i,
    //         //         jet.rap(),
    //         //         jet.phi(),
    //         //         jet.pt()
    //         //     )
    //         //     .unwrap();
    //         // }
    //         // let elapsed = now.elapsed()?.as_nanos();
    //         // times.push(elapsed);
    //         //}

    //         // let mean = times.iter().sum::<u128>() as f64 / times.len() as f64;

    //         // let variance_sum: f64 = times
    //         //     .iter()
    //         //     .map(|&x| {
    //         //         let diff = x as f64 - mean;
    //         //         diff * diff
    //         //     })
    //         //     .sum();

    //         // let sample_variance = variance_sum as f64 / (times.len() - 1) as f64;

    //         // let uncertainty = sample_variance.sqrt() / (times.len() as f64).sqrt();

    //         // writeln!(
    //         //     out,
    //         //     "The average time for 1000 iterations is {:.3} µs ± {:.3} µs",
    //         //     mean / 1000.0,
    //         //     uncertainty / 1000.0
    //         // )?;
    //     })
    // });
}

criterion_group!(benches, example_bench);
criterion_main!(benches);
