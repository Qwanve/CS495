use core::ops::ControlFlow;
use meshx::io::save_trimesh_ascii;
use rand::prelude::*;
use std::convert::identity;
use std::fmt::Display;
use std::fs::File;
use std::io::Write;
use std::num::NonZeroUsize;
use std::ops::{Add, Sub};
use std::time::{Duration, Instant};

use duration_string::DurationString;

use indicatif::{ProgressBar, ProgressStyle};

use clap::Parser;

use lazy_static::lazy_static;

use meshx::TriMesh;

// Parse the arguments and store them in a global variable
//  They never change, so a global is safe to use
lazy_static! {
    static ref ARGS: Args = Args::parse();
}

// Define what arguments the user can input.
#[derive(Parser, Debug)]
struct Args {
    /// Number of rings. Must be greater than 0
    #[arg()]
    rings: NonZeroUsize,

    /// Error margin for distance calculations
    #[arg(short, default_value_t = 10e-6)]
    error_margin: f64,

    // TODO: Better error reporting
    /// Time to run before exiting: [0-9]+(ns|us|ms|[smhdwy]). eg: 5m, 15s, 100m
    #[arg(short, default_value = "15s")]
    time: DurationString,

    /// Enable debug printing
    #[arg(short, action)]
    debug: bool,
}

impl Args {
    /// Convert the user given time string to a duration
    fn time(&self) -> Duration {
        self.time.into()
    }
}

// TODO: Change this to a vector class eg Vector3D
/// A Point type
#[derive(Copy, Clone)]
struct Point {
    x: f64,
    y: f64,
    z: f64,
    // Used for ring generation.
    is_cusp: bool,
}

impl Point {
    /// Create a new point. By default, points aren't cusps
    pub const fn new(x: f64, y: f64, z: f64) -> Point {
        Point {
            x,
            y,
            z,
            is_cusp: false,
        }
    }

    /// The distance between points
    /// # Example
    /// ```
    /// let point1 = Point::new(0.0, 1.0, 1.0);
    /// let point2 = Point::new(0.0, 0.0, 0.0);
    /// assert_eq!(point1.distance(&point2), 2.0.sqrt());
    ///```
    pub fn distance(&self, p: &Point) -> f64 {
        ((self.x - p.x).powi(2) + (self.y - p.y).powi(2) + (self.z - p.z).powi(2)).sqrt()
    }

    /// Scale a point as if it's a vector
    pub fn scale(&self, scalar: f64) -> Point {
        Point::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

// This allows adding two points together as if they're vectors
impl Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

// This allows subtracting two points as if they're vectors
impl Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

// Allows printing a point
impl Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{},{}", self.x, self.y, self.z)
    }
}

/// The base struct
/// This contains all the points, and keeps track of which points form lines
/// and triangles. Points are shared between triangles, so if you move one point,
/// the point in all triangles that share that point move as well
struct Mesh {
    points: Vec<Point>,
    pairs: Vec<(usize, usize)>,
    tris: Vec<(usize, usize, usize)>,
    ring_counts: Vec<usize>,
}

impl Mesh {
    /// Creates a new mesh with the proper amount of points
    pub fn new(rings: NonZeroUsize) -> Mesh {
        // Generate number of points per ring
        // Formula found https://oeis.org/A001354
        let mut fib = vec![0, 1];
        for i in 2..((usize::from(rings) + 1) * 2) {
            fib.push(fib[i - 1] + fib[i - 2]);
        }
        let mut ring_counts: Vec<usize> = fib.iter().step_by(2).map(|x| x * 7).collect();
        ring_counts[0] = 1;

        // Create random points
        let mut points: Vec<Point> = Vec::new();
        points.push(Point::new(0.0, 0.0, 0.0));
        for _ in 1..ring_counts.iter().sum() {
            points.push(sphere_rand(1.0));
        }

        let mut pairs: Vec<(usize, usize)> = Vec::new();
        let mut tris: Vec<(usize, usize, usize)> = Vec::new();

        // Manually prepare first ring to help the generation algorithm
        for i in 1..=7 {
            pairs.push((i, i % 7 + 1));
            pairs.push((0, i));
            tris.push((0, i, i % 7 + 1));
        }

        // generate every ring from 2 to n
        for ring in 2..=rings.into() {
            let offset: usize = ring_counts[..ring].iter().sum();
            let prev_offset: usize = ring_counts[..ring - 1].iter().sum();
            let mut cur_previous: usize = offset - 1;
            let mut to_next_cusp: usize = 0;
            for i in 0..ring_counts[ring] {
                let index = i + offset;
                if to_next_cusp == 0 {
                    points[index].is_cusp = true;
                    pairs.push((index, cur_previous));

                    let temp =
                        (cur_previous - prev_offset + 1) % ring_counts[ring - 1] + prev_offset;
                    pairs.push((index, temp));
                    tris.push((index, cur_previous, temp));
                    cur_previous = temp;
                    to_next_cusp = if points[cur_previous].is_cusp { 1 } else { 2 }
                } else {
                    pairs.push((index, cur_previous));
                    to_next_cusp -= 1;
                }
                let next_index = (index - offset + 1) % ring_counts[ring] + offset;
                pairs.push((index, next_index)); // ring lines
                tris.push((index, next_index, cur_previous));
            }
        }

        Mesh {
            points,
            pairs,
            tris,
            ring_counts,
        }
    }

    /// Does one iteration. For every point, it they are further apart
    ///  than 1, we move them closer together. If they are closer together
    ///  we push them further apart. If after checking every distance,
    ///  we move no points, then we return `ControlFlow::Break(())`.
    pub fn do_iteration(&mut self) -> ControlFlow<(), ()> {
        if self
            .pairs
            .iter()
            .copied()
            .map(|(a, b)| {
                if let ControlFlow::Continue((p1, p2)) =
                    move_points(&self.points[a], &self.points[b])
                {
                    self.points[a] = p1;
                    self.points[b] = p2;
                    false
                } else {
                    true
                }
            })
            //Bitwise-and to certainly prevent short-circuit behavior
            .fold(true, |acc, elem| acc & elem)
        {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    }

    /// Converts the mesh into a vector of points of the triangles
    pub fn get_tris(&self) -> Vec<(Point, Point, Point)> {
        self.tris
            .iter()
            .copied()
            .map(|(a, b, c)| (self.points[a], self.points[b], self.points[c]))
            .collect()
    }

    fn print_debug(&self) {
        for (a, b) in self.pairs.iter().copied() {
            println!(
                "Distance: {a}, {b} = {}",
                self.points[a].distance(&self.points[b])
            );
        }

        for (a, b, c) in self.tris.iter().copied() {
            let dist1 = self.points[a].distance(&self.points[b]);
            let dist2 = self.points[b].distance(&self.points[c]);
            let dist3 = self.points[c].distance(&self.points[a]);

            println!("Triangle: {a}, {b}, {c}");
            println!("Dists: {dist1}, {dist2}, {dist3}");
        }
    }
}

fn sphere_rand(radius: f64) -> Point {
    let mut rng = thread_rng();

    // Generate two random numbers between 0 and 1
    let u = rng.gen::<f64>();
    let v = rng.gen::<f64>();

    // Calculate the longitude and latitude
    let lon = 2.0 * std::f64::consts::PI * u;
    let lat = f64::acos(2.0_f64.mul_add(v, -1.0)) - std::f64::consts::FRAC_PI_2;

    // Calculate the x, y, and z coordinates of the point
    let x = f64::cos(lon) * f64::cos(lat) * radius;
    let y = f64::sin(lon) * f64::cos(lat) * radius;
    let z = f64::sin(lat) * radius;
    Point::new(x, y, z)
}

fn move_points(p1: &Point, p2: &Point) -> ControlFlow<(), (Point, Point)> {
    // Caculate the distance between the two points
    let dist = p1.distance(p2);
    // If the distance is close enough to the error margin
    if (dist - 1.0).abs() > ARGS.error_margin {
        return ControlFlow::Break(());
    }
    // Else, we move them by 1/10th the distance with a touch of randomness
    let scalar = (1.0 - dist) * 0.1;
    let offset1 = (*p1 - *p2).scale(scalar) + sphere_rand(scalar * 0.7);
    let offset2 = (*p2 - *p1).scale(scalar) + sphere_rand(scalar * 0.7);
    ControlFlow::Continue((*p1 + offset1, *p2 + offset2))
}

fn main() {
    println!("Creating mesh with {} rings", ARGS.rings);

    println!("Creating points");
    // Create the new mesh with the user given ring count
    let mut mesh = Mesh::new(ARGS.rings);

    println!("Moving points into proper place");
    // Start a timer to have time bounded ending
    let start = Instant::now();

    // Create a progress bar for the user
    let pb = ProgressBar::new(ARGS.time().as_millis() as u64);
    pb.set_style(
        ProgressStyle::with_template("{percent}%  {wide_bar}  [{elapsed_precise}]").unwrap(),
    );

    // While the timer hasn't surpassed the user given time-limit
    while start.elapsed() < ARGS.time() {
        // We move the points
        if let ControlFlow::Break(_) = mesh.do_iteration() {
            // If we didn't move any, we are within error margins and we can exit
            pb.println("Stopping as we are within error margins");
            break;
        }
        // Otherwise, update the progress bar, and continue
        pb.set_position(start.elapsed().as_millis() as u64);
    }
    pb.finish();

    // Debug print the mesh with every distance and point
    if ARGS.debug {
        mesh.print_debug();
    }

    // This was for outputting to csv, which we may come back to
    /*
    let mut output_file = File::create("./output.csv").unwrap();
    let tris = mesh.get_tris();
    for (a, b, c) in tris {
        let point1 = a.to_string();
        let point2 = b.to_string();
        let point3 = c.to_string();
        writeln!(output_file, "{point1:<40},{point2:<40},{point3:<40}").unwrap();
    }
    println!("Successfully created output.csv");
    */

    // Convert the Vec<Point> into Vec<[f64; 3]>
    let verts = mesh
        .points
        .iter()
        .copied()
        .map(|p| [p.x, p.y, p.z])
        .collect::<Vec<_>>();

    // Convert the Vec<(usize, usize, usize)> into Vec<[usize; 3]>
    let indices = mesh
        .tris
        .iter()
        .copied()
        .map(|(x, y, z)| [x, y, z])
        .collect::<Vec<_>>();

    // Create a mesh object
    let tri_mesh = TriMesh::new(verts, indices);
    // And save that object into an .obj file
    save_trimesh_ascii(&tri_mesh, "./output.obj").unwrap();
}
