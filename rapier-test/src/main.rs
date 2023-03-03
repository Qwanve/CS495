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

fn move_points(p1: &Point, p2: &Point) -> (Point, Point) {
    let dist = p1.distance(p2);
    if dist != 1.0 {
        let scalar = ((1.0 - dist)) * 0.1;
        let offset1 = (*p1 - *p2).scale(scalar) + sphere_rand(scalar * 0.7);
        let offset2 = (*p2 - *p1).scale(scalar) + sphere_rand(scalar * 0.7);
        return (*p1 + offset1, *p2 + offset2)
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
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();

    points.push(Point::new(0.0,0.0,0.0));

    for i in 1..=7 {
        points.push(sphere_rand(1.0));
        pairs.push((i, i%7+1));
        pairs.push((0, i));
        tris.push((0, i, i%7+1));
    }

    for i in 0..21 {
        let j = i+8;
        let j1 = (i+1)%21+8;
        points.push(sphere_rand(1.0));
        if (i%3 == 0) {
            pairs.push((j, (i/3+6)%7+1));
            pairs.push((j, (i/3)%7+1));
            tris.push((j, (i/3+6)%7+1, (i/3)%7+1));
        }
        else {
            pairs.push((j, i/3+1));
        }
        pairs.push((j, j1));
        tris.push((j, j1, i/3+1));
    }


    for _ in 0..200000 {
        for (a, b) in &pairs {
            let (p1, p2) = move_points(&points[*a], &points[*b]);
            points[*a] = p1;
            points[*b] = p2;
        }
    }

    for (a, b) in &pairs {
        println!("Distance: {}, {} = {}", *a, *b, points[*a].distance(&points[*b]));
    }

    for (a, b, c) in &tris {
        let dist1 = points[*a].distance(&points[*b]);
        let dist2 = points[*b].distance(&points[*c]);
        let dist3 = points[*c].distance(&points[*a]);

        println!("Triangle: {}, {}, {}", *a, *b, *c);
        println!("Dists: {}, {}, {}", dist1, dist2, dist3);
    }

}