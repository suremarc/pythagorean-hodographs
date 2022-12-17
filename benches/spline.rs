use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use glam::Vec3;
use pythagorean_hodographs::Spline;
use rand::{distributions::Standard, thread_rng, Rng};

fn spline(c: &mut Criterion) {
    c.bench_function("quintic PH spline (100 points)", |b| {
        b.iter_batched(
            || thread_rng().sample_iter(Standard).take(100).collect(),
            |data: Vec<Vec3>| black_box(Spline::catmull_rom(data)),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, spline);
criterion_main!(benches);
