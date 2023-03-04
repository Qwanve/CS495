use core::ops::ControlFlow;
use rand::prelude::*;
use std::convert::identity;
use std::fmt::Display;
use std::fs::File;
use std::io::Write;
use std::ops::{Add, Sub};
use std::time::{Duration, Instant};

use indicatif::{ProgressBar, ProgressStyle};

const ERROR_MARGIN: f32 = 10e-12;
const MAX_CALCULATION_TIME: Duration = Duration::from_secs(15);

#[derive(Copy, Clone)]
struct Point {
    x: f32,
    y: f32,
    z: f32,
}

impl Point {
    const fn new(x: f32, y: f32, z: f32) -> Point {
        Point { x, y, z }
    }

    fn distance(&self, p: &Point) -> f32 {
        ((self.x - p.x).powi(2) + (self.y - p.y).powi(2) + (self.z - p.z).powi(2)).sqrt()
    }

    fn scale(&self, scalar: f32) -> Point {
        Point::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{},{}", self.x, self.y, self.z)
    }
}

fn sphere_rand(radius: f32) -> Point {
    let mut rng = thread_rng();

    // Generate two random numbers between 0 and 1
    let u = rng.gen::<f32>();
    let v = rng.gen::<f32>();

    // Calculate the longitude and latitude
    let lon = 2.0 * std::f32::consts::PI * u;
    let lat = f32::acos(2.0_f32.mul_add(v, -1.0)) - std::f32::consts::FRAC_PI_2;

    // Calculate the x, y, and z coordinates of the point
    let x = f32::cos(lon) * f32::cos(lat) * radius;
    let y = f32::sin(lon) * f32::cos(lat) * radius;
    let z = f32::sin(lat) * radius;
    Point::new(x, y, z)
}

fn move_points(p1: &Point, p2: &Point) -> ControlFlow<(), (Point, Point)> {
    let dist = p1.distance(p2);
    if (dist - 1.0).abs() > ERROR_MARGIN {
        let scalar = (1.0 - dist) * 0.1;
        let offset1 = (*p1 - *p2).scale(scalar) + sphere_rand(scalar * 0.7);
        let offset2 = (*p2 - *p1).scale(scalar) + sphere_rand(scalar * 0.7);
        ControlFlow::Continue((*p1 + offset1, *p2 + offset2))
    } else {
        ControlFlow::Break(())
    }
}

fn main() {
    let mut points: Vec<Point> = Vec::new();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();

    println!("Creating points");
    points.push(Point::new(0.0, 0.0, 0.0));
    //Ring 1
    for i in 1..=7 {
        points.push(sphere_rand(1.0));
        pairs.push((i, i % 7 + 1));
        pairs.push((0, i));
        tris.push((0, i, i % 7 + 1));
    }

    for i in 0..21 {
        let j = i + 8;
        let j1 = (i + 1) % 21 + 8;
        points.push(sphere_rand(1.0));
        if i % 3 == 0 {
            pairs.push((j, (i / 3 + 6) % 7 + 1));
            pairs.push((j, (i / 3) % 7 + 1));
            tris.push((j, (i / 3 + 6) % 7 + 1, (i / 3) % 7 + 1));
        } else {
            pairs.push((j, i / 3 + 1));
        }
        pairs.push((j, j1));
        tris.push((j, j1, i / 3 + 1));
    }

    println!("Moving points into proper place");
    let start = Instant::now();
    let pb = ProgressBar::new(MAX_CALCULATION_TIME.as_millis() as u64);
    pb.set_style(
        ProgressStyle::with_template("{percent}%  {wide_bar}  [{elapsed_precise}]").unwrap(),
    );
    while start.elapsed() < MAX_CALCULATION_TIME {
        let count = pairs
            .iter()
            .copied()
            .map(|(a, b)| {
                if let ControlFlow::Continue((p1, p2)) = move_points(&points[a], &points[b]) {
                    points[a] = p1;
                    points[b] = p2;
                    false
                } else {
                    true
                }
            })
            .collect::<Vec<bool>>();
        if count.into_iter().all(identity) {
            println!("Stopping as we are within error margins");
            break;
        }
        pb.set_position(start.elapsed().as_millis() as u64);
    }
    pb.finish();

    for (a, b) in pairs.iter().copied() {
        println!("Distance: {a}, {b} = {}", points[a].distance(&points[b]));
    }

    for (a, b, c) in tris.iter().copied() {
        let dist1 = points[a].distance(&points[b]);
        let dist2 = points[b].distance(&points[c]);
        let dist3 = points[c].distance(&points[a]);

        println!("Triangle: {a}, {b}, {c}");
        println!("Dists: {dist1}, {dist2}, {dist3}");
    }
    let mut output_file = File::create("./output.csv").unwrap();
    for (a, b, c) in tris.iter().copied() {
        let point1 = points[a].to_string();
        let point2 = points[b].to_string();
        let point3 = points[c].to_string();
        writeln!(output_file, "{point1:<40},{point2:<40},{point3:<40}").unwrap();
    }
    println!("Successfully created output.csv");
}
