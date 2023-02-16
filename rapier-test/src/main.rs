use rand::prelude::*;
use std::ops::{Add, Sub};

#[derive(Copy, Clone)]
struct Point {
    x: f32,
    y: f32,
    z: f32
}

impl Point {
    fn new(x: f32, y: f32, z: f32) -> Point {
        Point{x, y, z}
    }

    fn distance(&self, p: &Point) -> f32 {
        ((self.x - p.x).powf(2.0) + (self.y - p.y).powf(2.0) + (self.z - p.z).powf(2.0)).sqrt()
    }

    fn scale(&self, scalar: f32) -> Point{
        Point::new(self.x*scalar, self.y*scalar, self.z*scalar)
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
        Point::new(self.x+other.x, self.y+other.y, self.z+other.z)
    }
}

impl Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point::new(self.x-other.x, self.y-other.y, self.z-other.z)
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

fn main() {
    let center = Point::new(0.0,0.0,0.0);

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
    }

    for i in 0..7 {
        let i2 = (i+1)%7;
        println!("Point {}: {}, {}, {}", i, ring1[i].x, ring1[i].y, ring1[i].z);
        println!("Distance to next: {}", ring1[i].distance(&ring1[i2]));
        println!("Distance to center: {}\n", ring1[i].distance(&center));
    }
}