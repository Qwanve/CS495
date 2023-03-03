use rand::prelude::*;
use std::{
    fmt::Display,
    ops::{Add, Sub},
};

//Should this derive partialeq?
#[derive(Copy, Clone, PartialEq)]
struct Point {
    x: f32,
    y: f32,
    z: f32,
}

impl Point {
    fn new(x: f32, y: f32, z: f32) -> Point {
        Point { x, y, z }
    }

    fn distance(&self, p: &Point) -> f32 {
        ((self.x - p.x).powf(2.0) + (self.y - p.y).powf(2.0) + (self.z - p.z).powf(2.0)).sqrt()
    }

    fn scale(&self, scalar: f32) -> Point {
        Point::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}, {}, {}", self.x, self.y, self.z)
    }
}

impl From<(f32, f32, f32)> for Point {
    fn from(t: (f32, f32, f32)) -> Point {
        let (x, y, z) = t;
        Point::new(x, y, z)
    }
}

impl From<Point> for (f32, f32, f32) {
    fn from(p: Point) -> (f32, f32, f32) {
        (p.x, p.y, p.z)
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

fn sphere_rand(radius: f32) -> Point {
    let mut rng = thread_rng();

    // Generate two random numbers between 0 and 1
    let u = rng.gen::<f32>();
    let v = rng.gen::<f32>();

    // Calculate the longitude and latitude
    let lon = 2.0 * std::f32::consts::PI * u;
    let lat = f32::acos(2.0 * v - 1.0) - std::f32::consts::FRAC_PI_2;

    // Calculate the x, y, and z coordinates of the point
    let x = f32::cos(lon) * f32::cos(lat) * radius;
    let y = f32::sin(lon) * f32::cos(lat) * radius;
    let z = f32::sin(lat) * radius;
    Point::new(x, y, z)
}

fn move_points(p1: &Point, p2: &Point) -> (Point, Point) {
    let dist = p1.distance(p2);
    if dist != 1.0 {
        let scalar = (1.0 - dist) * 0.1;
        let offset1 = (*p1 - *p2).scale(scalar) + sphere_rand(scalar * 0.7);
        let offset2 = (*p2 - *p1).scale(scalar) + sphere_rand(scalar * 0.7);
        return (*p1 + offset1, *p2 + offset2);
    }
    (*p1, *p2)
}

fn main() {
    /*let center = Point::new(0.0,0.0,0.0);

    let mut ring1: Vec<Point> = Vec::new();

    for _ in 0..7 {
        ring1.push(sphere_rand(1.0));
    }

    for _ in 0..100000 {
        for i in 0..7 {
            let i2 = (i+1)%7;
            let distc = ring1[i].distance(&center);
            let offset = (ring1[i].scale(1.0/distc) - ring1[i]).scale(0.5);
            ring1[i] = ring1[i] + offset;

            let dist = ring1[i].distance(&ring1[i2]);
            //println!("{}", dist);
            if dist != 1.0 {
                let scalar = ((1.0 - dist)) * 0.1;
                let offset1 = (ring1[i] - ring1[i2]).scale(scalar) + sphere_rand(scalar * 0.7);
                let offset2 = (ring1[i2] - ring1[i]).scale(scalar) + sphere_rand(scalar * 0.7);
                ring1[i] = ring1[i] + offset1;
                ring1[i2] = ring1[i2] + offset2;
            }
        }
    }*/

    let mut points: Vec<Point> = Vec::new();
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    points.push(Point::new(0.0, 0.0, 0.0));

    for i in 1..=7 {
        points.push(sphere_rand(1.0));
        pairs.push((i, (i + 1) % 7));
        pairs.push((0, i));
    }

    for i in 0..21 {
        let j = i + 8;
        let j1 = (i + 1) % 21 + 8;
        points.push(sphere_rand(1.0));
        if i % 3 == 0 {
            pairs.push((j, i / 3));
            pairs.push((j, (i / 3 + 1) % 7));
        } else {
            pairs.push((j, i / 3));
        }
        pairs.push((j, j1));
    }

    for _ in 0..100000 {
        for (a, b) in &pairs {
            let (p1, p2) = move_points(&points[*a], &points[*b]);
            points[*a] = p1;
            points[*b] = p2;
        }
    }
    println!("Finished calculating points");
    println!("Finding triangles");

    // for (a, b) in &pairs {
    //     println!("Distance: {}", points[*a].distance(&points[*b]));
    // }

    let mut triangles = vec![];
    for (a, b) in pairs.iter().copied() {
        for (c, d) in pairs.iter().copied() {
            if (a, b) == (c, d) {
                continue;
            }
            let (common, third) = match (a, b, c, d) {
                _ if a == c => (a, (b, d)),
                _ if a == d => (a, (b, c)),
                _ if b == c => (b, (a, d)),
                _ if b == d => (b, (a, c)),
                _ => continue,
            };
            if triangles.iter().any(|t: &Vec<Point>| {
                t.iter()
                    .all(|&p| p == points[common] || p == points[third.0] || p == points[third.1])
            }) {
                continue;
            }
            if pairs
                .iter()
                .find(|&&x| x == third || x == (third.1, third.0))
                .is_some()
            {
                let point_1 = points[common];
                let point_2 = points[third.0];
                let point_3 = points[third.1];
                let triangle = vec![point_1, point_2, point_3];
                triangles.push(triangle);
            }
        }
    }
    println!("Found {} triangles", triangles.len());
    for t in triangles {
        println!("Triangle:  {},  {},  {}", t[0], t[1], t[2]);
    }
}
