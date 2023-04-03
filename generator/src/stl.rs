use stl_io::{Vector, Normal, Triangle};
use super::Point;

impl From<Point> for stl_io::Vertex {
    fn from(value: Point) -> Self {
        Vector::new([value.x as f32, value.y as f32, value.z as f32])
    }
}

fn normal((x, y, z): (Point, Point, Point)) -> Normal {
    let leg_a = x - y;
    let leg_b = z - x;
    let normal_x = leg_a.y * leg_b.z - leg_a.z * leg_b.y;
    let normal_y = leg_a.z * leg_b.x - leg_a.x * leg_b.z;
    let normal_z = leg_a.x * leg_b.y - leg_a.y * leg_b.x;
    Normal::new([normal_x as f32, normal_y as f32, normal_z as f32])
}

pub(crate) fn to_stl_triangle((x, y, z): (Point, Point, Point)) -> Triangle {
    let normal = normal((x, y, z));
    Triangle {
        normal,
        vertices: [
            x.into(),
            y.into(),
            z.into()
        ]
    }
}