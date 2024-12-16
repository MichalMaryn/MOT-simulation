extern crate atomecs as lib;
extern crate nalgebra;
use lib::atom::{Atom, Position, Velocity};
use lib::atom_sources::emit::AtomNumberToEmit;
use lib::atom_sources::oven::{OvenAperture, OvenBuilder};
use lib::atom_sources::mass::{MassDistribution, MassRatio};
use lib::atom_sources::{VelocityCap, AtomSourcePlugin};
use lib::destructor::ToBeDestroyed;
use lib::integrator::Timestep;
use lib::laser::LaserPlugin;
use lib::laser::gaussian::GaussianBeam;
use lib::laser_cooling::{CoolingLight, LaserCoolingPlugin};
use lib::magnetic::quadrupole::QuadrupoleField3D;
use lib::output::file::FileOutputPlugin;
use lib::output::file::Text;
use lib::shapes::{Cuboid, Cylinder};
use lib::gravity::ApplyGravitationalForceSystem;
use lib::gravity::ApplyGravityOption;
use lib::sim_region::{SimulationVolume, VolumeType};
use lib::simulation::SimulationBuilder;
use lib::species::{Silver109, Silver109_328};
use nalgebra::Vector3;
use specs::prelude::*;
use std::time::Instant;
extern crate specs;
#[allow(unused_imports)]
use specs::{Builder, Entity, RunNow, World};
#[allow(unused_imports)]

const BEAM_NUMBER : usize = 13;

pub fn calculate_rayleigh_range(wavelength: &f64, e_radius: &f64) -> f64 {
    2.0 * 3.14 * e_radius.powf(2.0) / wavelength
}

fn main() {
    let now = Instant::now();

    let mut sim_builder = SimulationBuilder::default();
    sim_builder.add_plugin(LaserPlugin::<{BEAM_NUMBER}>);
    sim_builder.add_plugin(LaserCoolingPlugin::<Silver109_328, {BEAM_NUMBER}>::default());

    sim_builder.add_plugin(AtomSourcePlugin::<Silver109>::default());
    sim_builder.add_plugin(FileOutputPlugin::<Position, Text, Atom>::new("pos.txt".to_string(), 100));
    sim_builder.add_plugin(FileOutputPlugin::<Velocity, Text, Atom>::new("vel.txt".to_string(), 100));
    let mut sim = sim_builder.build();

    // Magnetic gradient
    let quadrupole_gradient_2d = 12.885;                    
    let quadrupole_gradient_3d = 8.5;
                                             
    // Create cooling lasers - push beam.
    let detuning_push = -63.72;
    let power_push = 0.49234;
    let radius_push = 1.6178;

    // Create cooling lasers - 2D MOT.
    let detuning_2d = detuning_push;
    let power_2d = 249.0;                                                      
    let radius_2d = 12.442;

    let ellipticity = 0.56;

    // Create cooling lasers - 3D MOT.
    let detuning_3d = detuning_push;
    let power_3d = 50.0;
    let radius_3d = 7.5;

    let size_x_oven = 18.278;
    let size_y_oven = 1.0;

    let shift_z = 2.0; // 2.5
    let distance_oven_2d = 15.0; // minimal distance_oven_2d is 10.68 mm with radius_2d 8 mm and ellipticity 0.667
    let distance_2d_3d = 210.0;

    let atom_number = 10000000;
    let temperature = 900.0 + 273.15;
    let velocity_cap = 55.0;
    let time = 80;

    let dir_z = Vector3::new(0.0, 0.0, 1.0).normalize();
    let perp_x_z = dir_z.normalize().cross(&dir_z);
    let perp_y_z = dir_z.normalize().cross(&perp_x_z);

    let dir_x = Vector3::new(1.0, 0.0, 0.0).normalize();
    let perp_x_x = dir_x.normalize().cross(&dir_x);
    let perp_y_x = dir_x.normalize().cross(&perp_x_x);

    // Best observed feasible point: 04/06/2024 best 135/5000000
//     Best estimated feasible point (according to models):
//     ellipticity    detuning    ratio_power_2d_3d    radius_2d    quadrupole_gradient_2d    quadrupole_gradient_3d
//     ___________    ________    _________________    _________    ______________________    ______________________

//       56.151        -63.72          83.289           12.442              12.885                    5.1218        

// Estimated objective function value = -3.058e-05

    sim.world.insert(ApplyGravityOption);

    //Create magnetic field 2D MOT.
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(quadrupole_gradient_2d, Vector3::x()))
        .with(Position {pos: Vector3::new(0.0, 0.0, 0.0)})
        .build();

    //Create magnetic field 3D MOT.
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(quadrupole_gradient_3d, Vector3::z()))
        .with(Position {pos: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001)})
        .build();

    // 1 push beam
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_push * 0.001 / 2.0_f64.sqrt(),
            power: power_push * 0.001,
            direction: Vector3::new(1.0, 0.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_push,
            -1,
        ))
        .build();
    
    // 4 cooling lasers of 2D MOT
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 1.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, -1.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 1.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, -1.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_2d,
            1,
        ))
        .build();

    // 6 cooling lasers of 3D MOT
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(1.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(-1.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(-1.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(1.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(0.0, 0.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            -1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance_2d_3d * 0.001, 0.0, -shift_z * 0.001),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001 / 3.0,
            direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            -1,
        ))
        .build();

    // Create an oven - effusion cell
    // The oven will eject atoms on the first frame and then be deleted.
    // Oven representing the dispenser

    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Silver109>::new(
                temperature,
                Vector3::new(0.0, 0.0, 1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [size_x_oven * 0.001,size_y_oven * 0.001, 0.000127]}) //0.000564189
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, 0.0, -(distance_oven_2d * 0.001 + 0.000127/2.0)),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 109.0,
            ratio: 1.0,
        }]))
        .with(AtomNumberToEmit {
            number: atom_number as i32,
        })
        .with(ToBeDestroyed)
        .build();

    // Use a simulation bound so that atoms that escape the capture region are deleted from the simulation

    // Oven
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, -(distance_oven_2d * 0.001 + 0.000127/2.0)),
        })
        .with(Cuboid {
            half_width: Vector3::new(size_x_oven * 0.001 / 2.0, size_y_oven * 0.001 / 2.0, 0.000127 / 2.0),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //Tubed aparature
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, -((distance_oven_2d * 0.001)/2.0)),
        })
        .with(Cylinder {
            radius: 16.0e-3,
            length: distance_oven_2d * 0.001,
            direction: dir_z,
            perp_x: perp_x_z,
            perp_y: perp_y_z
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // 2d MOT chamber
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.02, 0.02, 0.02),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Small aperture leading through the mirror
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.02 + 0.01/2.0, 0.0, 0.0),
        })
        .with(Cylinder {
            radius: 0.00075,
            length: 0.01,
            direction: dir_x,
            perp_x: perp_x_x,
            perp_y: perp_y_x
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Wider tube leading to 3D MOT
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.02 + 0.01 + (distance_2d_3d * 0.001 - 0.01 - 0.01 - 0.02)/2.0, 0.0, 0.0),
        })
        .with(Cylinder {
            radius: 0.0035,
            length: (distance_2d_3d * 0.001 - 0.01 - 0.01 - 0.02),
            direction: dir_x,
            perp_x: perp_x_x,
            perp_y: perp_y_x
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // 3d MOT chamber
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(distance_2d_3d * 0.001, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.01, 0.01, 0.01),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Also use a velocity cap so that fast atoms are not even simulated.
    sim.world.insert(VelocityCap { value: velocity_cap });

    // Define timestep
    sim.world.insert(Timestep { delta: 1.0e-6 });

    sim.world.insert(ApplyGravitationalForceSystem);

    // Run the simulation for a number of steps.
    for _i in 0..(time * 1000) {
        sim.step();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}