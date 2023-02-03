use plotpy::{Canvas, Plot};
fn main() -> Result<(), &'static str> {
    let mut canvas = Canvas::new();
    add_triangle(
        &mut canvas,
        [0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [5.0, 0.0, 0.0],
    );
    let mut plot = Plot::new();
    plot.add(&canvas);
    plot.set_range_3d(0.0, 15.0, 0.0, 15.0, 0.0, 15.0)
        .set_frame_borders(false)
        .set_hide_axes(true)
        .set_equal_axes(true)
        .set_show_errors(true);
    plot.save_and_show("test.svg")?;
    Ok(())
}

fn add_triangle(c: &mut Canvas, point1: [f32; 3], point2: [f32; 3], point3: [f32; 3]) {
    c.draw_polyline(&[point1, point2, point3], true);
}
