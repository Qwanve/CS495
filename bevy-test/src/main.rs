use bevy::pbr::wireframe::{Wireframe, WireframePlugin};
use bevy::prelude::*;
use bevy::render::mesh::{self, PrimitiveTopology};

use serde::Deserialize;
use serde::Serialize;

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(WireframePlugin)
        .add_startup_system(setup)
        .run();
}

type Point = [f32; 3];
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
struct Triangle(Point, Point, Point);

impl Triangle {
    fn points(self) -> [[Point; 3]; 2] {
        [[self.0, self.1, self.2], [self.0, self.2, self.1]]
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
		.comment(Some(b'#'))
        .trim(csv::Trim::All)
        .from_path("./points.csv")
        .unwrap();
    let verts = rdr
        .deserialize::<Triangle>()
        .map(Result::unwrap)
        .map(Triangle::points)
        .collect::<Vec<[[[f32; 3]; 3]; 2]>>()
        .concat() //Vec<[[f32; 3]; 3]>
        .concat(); //Vec<[f32; 3]>

    let indices: Vec<u32> = (0..verts.len() as u32).collect();

    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, verts);

    mesh.set_indices(Some(mesh::Indices::U32(indices)));

    commands.spawn((
        PbrBundle {
            mesh: meshes.add(mesh),
            material: materials.add(Color::rgba(0.3, 0.5, 0.3, 0.0).into()),
            ..default()
        },
        Wireframe,
    ));

    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 1500.0,
            shadows_enabled: true,
            ..default()
        },
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(-2.0, -2.0, -2.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}
