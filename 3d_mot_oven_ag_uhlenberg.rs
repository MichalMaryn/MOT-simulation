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

const BEAM_NUMBER : usize = 8;

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

    sim.world.insert(ApplyGravityOption);

    // Shift in z axis between the 3d and 3D MOT and the x axis distance between them 'It equals 347 mm and includes a 55.7 mm-long differential pumping section, which is inserted into the 3d-MOT chamber'
    let shift = 0.0;

    // Magnetic gradient
    let magnetic_grad_3d = 10.0;

    // // Create cooling lasers - 3dMOT.
    let detuning_3d = -1.0 * 23.4;
    let power_3d = 1.5;
    let radius_3d = 3.0;

    let distance = 500.0;

    let number_to_emit = 2147483647;
    let temperature = 300.00;

    let vel_cap = 40.0;
    let time = 80;

    let dir_y = Vector3::new(0.0, 1.0, 0.0).normalize();
    let perp_x_y = dir_y.normalize().cross(&dir_y);
    let perp_y_y = dir_y.normalize().cross(&perp_x_y);

    //Create magnetic field 3D MOT.
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(magnetic_grad_3d, Vector3::z()))
        .with(Position {pos: Vector3::new(0.0, 0.0, 0.0)})
        .build();

    // 6 cooling lasers of 3D MOT
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
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
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
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
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
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
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
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
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
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
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_3d * 0.001,
            direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Silver109_328>(
            detuning_3d,
            -1,
        ))
        .build();

    // Add a 7th laser to fix a bug with a limiting oven
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, 0.0),
            e_radius: 0.5,
            power: 1.0e-30,
            direction: Vector3::new(0.0, -1.0, 0.0).normalize(),
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
        .with(Position {
            pos: Vector3::new(0.0, -(distance * 0.001 + 0.001/2.0), shift),
        })
        .with(Cylinder {
            radius: 0.000564189,
            length: 0.001,
            direction: dir_y,
            perp_x: perp_x_y,
            perp_y: perp_y_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Silver109>::new(
                temperature,
                Vector3::new(0.0, 1.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [0.000564189*2.0,0.001,0.000564189*2.0]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, -(distance * 0.001 + 0.001/2.0), shift),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 109.0,
            ratio: 1.0,
        }]))
        .with(AtomNumberToEmit {
            number: number_to_emit,
        })
        .with(ToBeDestroyed)
        .build();

    // Use a simulation bound so that atoms that escape the capture region are deleted from the simulation
    // Tube from oven to hot lip. - tubed aparature
    //Tubed aparature
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, -((distance * 0.001 - 0.006)/2.0 + 0.006), shift),
        })
        .with(Cylinder {
            radius: 10.0e-3,
            length: 0.494,
            direction: dir_y,
            perp_x: perp_x_y,
            perp_y: perp_y_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Endcap hole leading to the chamber
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, -(0.005 + 0.0005), shift),
        })
        .with(Cylinder {
            radius: 10.0e-3,
            length: 0.001,
            direction: dir_y,
            perp_x: perp_x_y,
            perp_y: perp_y_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //Aparature filtering atoms that will have interaction with atoms
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, -(0.0045 + 0.0005/2.0), shift),
        })
        .with(Cylinder {
            radius: 3.4e-3,
            length: 0.0005,
            direction: dir_y,
            perp_x: perp_x_y,
            perp_y: perp_y_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // 3d MOT chamber
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.0045, 0.0045, 0.0045),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Also use a velocity cap so that fast atoms are not even simulated.
    sim.world.insert(VelocityCap { value: vel_cap });

    // Define timestep
    sim.world.insert(Timestep { delta: 1.0e-6 });

    sim.world.insert(ApplyGravitationalForceSystem);

    // Run the simulation for a number of steps.
    for _i in 0..(time * 1000) {
        sim.step();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}