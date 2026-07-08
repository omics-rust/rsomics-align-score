use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-align-score"))
}

#[test]
fn identical_sequences_identity_1() {
    let dir = std::env::temp_dir().join("align-compat-ident");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("ident.fa");
    std::fs::write(&fa, ">a\nATCGATCG\n>b\nATCGATCG\n").unwrap();

    let out = Command::new(ours()).arg(&fa).output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let line = s.lines().nth(1).unwrap();
    let score: i32 = line.split('\t').nth(2).unwrap().parse().unwrap();
    let identity: f64 = line.split('\t').nth(3).unwrap().parse().unwrap();

    assert_eq!(score, 8, "8bp identical with match=1 should score 8");
    assert!(
        (identity - 1.0).abs() < 0.001,
        "identical seqs should have identity=1.0, got {identity}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn score_increases_with_similarity() {
    let dir = std::env::temp_dir().join("align-compat-sim");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let similar = dir.join("sim.fa");
    std::fs::write(&similar, ">a\nATCGATCG\n>b\nATCGATCG\n").unwrap();
    let different = dir.join("diff.fa");
    std::fs::write(&different, ">a\nATCGATCG\n>b\nTTTTTTTT\n").unwrap();

    let out_sim = Command::new(ours()).arg(&similar).output().unwrap();
    let out_diff = Command::new(ours()).arg(&different).output().unwrap();

    let score_sim: i32 = String::from_utf8(out_sim.stdout)
        .unwrap()
        .lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .parse()
        .unwrap();
    let score_diff: i32 = String::from_utf8(out_diff.stdout)
        .unwrap()
        .lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .parse()
        .unwrap();

    assert!(
        score_sim > score_diff,
        "identical should score higher ({score_sim}) than different ({score_diff})"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// Value-exact to Bio.Align.PairwiseAligner (match=1, mismatch=-1,
/// open_gap_score=-2, extend_gap_score=-1) on gapped alignments — the case the
/// original ungapped-only tests never covered. Golden scores from Biopython
/// 1.87. A gap of length k costs -2-(k-1); scores were verified over hundreds
/// of random pairs against Bio.Align.
#[test]
fn gapped_scores_match_bio_align() {
    let dir = std::env::temp_dir().join("align-compat-gapped");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let score = |file: &std::path::Path, local: bool| -> (i32, f64) {
        let mut cmd = Command::new(ours());
        cmd.arg(file);
        if local {
            cmd.arg("--local");
        }
        let out = cmd.output().unwrap();
        assert!(out.status.success());
        let s = String::from_utf8(out.stdout).unwrap();
        let row = s.lines().nth(1).unwrap();
        let f: Vec<&str> = row.split('\t').collect();
        (f[2].parse().unwrap(), f[3].parse().unwrap())
    };

    // 4-base leading gap in b: global keeps the gap (-5), local trims it.
    let lead = dir.join("lead.fa");
    std::fs::write(&lead, ">a\nGGGGACGTACGT\n>b\nACGTACGT\n").unwrap();
    assert_eq!(score(&lead, false).0, 3, "global leading-gap score");
    assert_eq!(score(&lead, true).0, 8, "local should trim the leading gap");

    // single internal insertion in b.
    let ins = dir.join("ins.fa");
    std::fs::write(&ins, ">a\nACGTACGTAC\n>b\nACGTTACGTAC\n").unwrap();
    let (g, ident) = score(&ins, false);
    assert_eq!(g, 8, "global single-insertion score");
    assert!((ident - 10.0 / 11.0).abs() < 0.001, "identity {ident}");

    let _ = std::fs::remove_dir_all(&dir);
}
