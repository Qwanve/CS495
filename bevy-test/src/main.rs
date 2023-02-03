use std::f32::consts::PI;
use bevy::prelude::*;
use bevy::render::mesh::{self, PrimitiveTopology};
use bevy::pbr::wireframe::{Wireframe, WireframeConfig, WireframePlugin};
use bevy::render::{render_resource::WgpuFeatures, settings::WgpuSettings, RenderPlugin};

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(WireframePlugin)
        .add_startup_system(setup)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList);

    // Positions of the vertices
    // See https://bevy-cheatbook.github.io/features/coords.html
    /*mesh.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        vec![[0., 0., 0.], [-0.5, 0.866, 0.], [0.5, 0.866, 0.]],
    );// */

    let verts = vec![[0., 0., 0.], [-0.5, 0.866, 0.], [-1., 0., 0.], [0., 0., 0.], [0.5, 0.866, 0.], [-0.5, 0.866, 0.]];
    let indices: Vec<u32> = (0..verts.len() as u32).collect();

    mesh.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        verts,
    );// */

    // In this example, normals and UVs don't matter,
    // so we just use the same value for all of them
    //mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, vec![[0., 1., 0.]; 6]);
    //mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, vec![[0., 0.]; 6]);

    // A triangle using vertices 0, 2, and 1.
    // Note: order matters. [0, 1, 2] will be flipped upside down, and you won't see it from behind!
    //mesh.set_indices(Some(mesh::Indices::U32(vec![0,1,2,3,4,5])));
    mesh.set_indices(Some(mesh::Indices::U32(indices)));

    commands.spawn((PbrBundle {
        mesh: meshes.add(mesh),
        material: materials.add(Color::rgba(0.3, 0.5, 0.3, 0.0).into()),
        ..default()
    }, Wireframe, ));

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
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}