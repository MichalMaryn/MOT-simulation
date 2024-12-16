// A 2D+ mot configuration, loaded directly from oven then a 3D MOT

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
use lib::laser::gaussian::CircularMask;
use lib::laser_cooling::{CoolingLight, LaserCoolingPlugin};
// use lib::laser_cooling::force::{EmissionForceConfiguration, EmissionForceOption};
// use lib::laser_cooling::photons_scattered::ScatteringFluctuationsOption;
use lib::magnetic::quadrupole::QuadrupoleField3D;
use lib::output::file::FileOutputPlugin;
use lib::output::file::Text;
use lib::shapes::{Cuboid, Cylinder};
use lib::gravity::ApplyGravitationalForceSystem;
use lib::gravity::ApplyGravityOption;
use lib::sim_region::{SimulationVolume, VolumeType};
use lib::simulation::SimulationBuilder;
use lib::species::{Cesium133, Cesium133_852D2};
use nalgebra::Vector3;
use specs::prelude::*;
use std::time::Instant;
extern crate specs;
#[allow(unused_imports)]
use specs::{Builder, Entity, RunNow, World};
#[allow(unused_imports)]

const BEAM_NUMBER : usize = 12;

pub fn calculate_rayleigh_range(wavelength: &f64, e_radius: &f64) -> f64 {
    2.0 * 3.14 * e_radius.powf(2.0) / wavelength
}

fn main() {
    let now = Instant::now();

    let mut sim_builder = SimulationBuilder::default();
    sim_builder.add_plugin(LaserPlugin::<{BEAM_NUMBER}>);
    sim_builder.add_plugin(LaserCoolingPlugin::<Cesium133_852D2, {BEAM_NUMBER}>::default());
    sim_builder.add_plugin(AtomSourcePlugin::<Cesium133>::default());
    sim_builder.add_plugin(FileOutputPlugin::<Position, Text, Atom>::new("pos.txt".to_string(), 100));
    sim_builder.add_plugin(FileOutputPlugin::<Velocity, Text, Atom>::new("vel.txt".to_string(), 100));
    // sim_builder.add_plugin(FileOutputPlugin::<MagneticFieldSampler, Text, Atom>::new("mag.txt".to_string(), 1));
    let mut sim = sim_builder.build();

    sim.world.insert(ApplyGravityOption);
// 2mm: 284.55        0.75524        -25.754         17.767 2.5mm: 299	      0.69      	 -25	         18.87 3mm: 241.64        0.28328        -22.994         19.692 4mm: 298.43        0.799          -25.239         17.665 5mm: 299.22        0.57828        -25.386         13.191 6mm 285.11         0.70926        -26.056         8.552 7.5mm: 285         0.65289        -25.218         5.0792

// Current best when B_3D = 11.05
// Trap rate in 3D MOT: 1.38696697895938*10^7 
// Flux from 2D MOT: 9.04729229351968*10^7 
// FEfficiency of 2D MOT: 0.1625*10^(-5) 
// When B=13
// Trap rate in 3D MOT: 0.61880065215111*10^7 
// Flux from 2D MOT: 8.68454708363799*10^7 
// FEfficiency of 2D MOT: 0.0725*10^(-5) 

// Magnetic gradient
    let quadrupole_gradient_2d = 13.0;
    let quadrupole_gradient_3d = 13.0;

    // Push beam along x axis
    let detuning_push= -1.9*5.22; //-26.232
    let power_push = 3.0*0.75;
    let radius_push = 2.2; // 1.5mm 1/e^2 radius because the hole limits it
    let wavelenght= 852.0e-9;

    // // Create cooling lasers - 2DMOT
    let detuning_2d = detuning_push;
    let power_2d = 360.0;
    let radius_2d = 34.0;
    let ellipticity = 0.972;

    let detuning_3d = detuning_push;
    // let power_3d = 95.0;
    let radius_3d = 6.9;

    let ratio_power_3d_2d = 0.8333;

    let number_to_emit = 1000000000.0;
    let velocity_cap = 70.0;
    let temperature = 300.0;

    let distance = 610.0;
    let shift_z = 3.0;
    let thick = 1.0;

    let time = 100;
    let small_radius = ((-1.0 * (ellipticity * ellipticity * radius_2d * 0.001 * radius_2d * 0.001 - radius_2d * 0.001 * radius_2d * 0.001)) as f64).sqrt();

    let number_to_emit_b = (number_to_emit*(small_radius*2.0*small_radius*2.0/(2.0*small_radius*2.0*small_radius*2.0 + 4.0*small_radius*2.0*radius_2d*0.001*2.0))) as i32;
    let number_to_emit_s = (number_to_emit*(small_radius*2.0*radius_2d*0.001*2.0/(2.0*small_radius*2.0*small_radius*2.0 + 4.0*small_radius*2.0*radius_2d*0.001*2.0))) as i32;

    let dir = Vector3::new(1.0, 0.0, 0.0).normalize();
    let perp_x = dir.normalize().cross(&dir);
    let perp_y = dir.normalize().cross(&perp_x);

    // Adding the 2D MOT magnetic field in x axis
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(quadrupole_gradient_2d, Vector3::x()))
        .with(Position {pos: Vector3::new(0.0, 0.0, shift_z * 0.001)})
        .build();

    //Create magnetic field 3D MOT.
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(quadrupole_gradient_3d, Vector3::z()))
        .with(Position {pos: Vector3::new(distance * 0.001, 0.0, 0.0)})
        .build();

    // Add push beam
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift_z * 0.001),
            e_radius: radius_push * 0.001 / 2.0_f64.sqrt(),
            power: power_push * 0.001,
            direction: Vector3::new(1.0,0.0,0.0).normalize(),
            rayleigh_range: calculate_rayleigh_range(&wavelenght, &(radius_push * 0.001 / 2.0_f64.sqrt())),
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_push,
            -1,
        ))
        .with(CircularMask{
            radius: 0.00075
        })
        .build();

    // Add 4 2D MOT lasers
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 0.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_2d,
            1,
        ))
        .build();

    // 6 cooling lasers of 3D MOT
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 6.0,
            direction: Vector3::new(1.0, 0.41, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 6.0,
            direction: Vector3::new(-1.0, -0.41, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 6.0,
            direction: Vector3::new(-0.41, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 6.0,
            direction: Vector3::new(0.41, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 3.0,
            direction: Vector3::new(0.0, 0.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            -1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(distance * 0.001, 0.0, 0.0),
            e_radius: radius_3d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 3.0,
            direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Cesium133_852D2>(
            detuning_3d,
            -1,
        ))
        .build();

    // Create an oven.
    // The oven will eject atoms on the first frame and then be deleted.
    // Atoms are taken from background gas which is inserted by 6 'ovens' - upper big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(0.0, 0.0, -1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, small_radius*2.0, thick * 0.001]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, 0.0, small_radius),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - lower big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(0.0, 0.0, 1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, small_radius*2.0, thick * 0.001]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, 0.0, -small_radius),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - back big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(0.0, -1.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, thick * 0.001, small_radius*2.0]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, small_radius, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - front big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(0.0, 1.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, thick * 0.001, small_radius*2.0]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, -small_radius, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - left small plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(1.0, 0.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [thick * 0.001, small_radius * 2.0, small_radius * 2.0]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(-radius_2d * 0.001, 0.0, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - right small plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Cesium133>::new(
                temperature,
                Vector3::new(-1.0, 0.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [thick * 0.001, small_radius*2.0, small_radius*2.0]})
                .with_microchannels(1.0, 0.0015) 
                .build(),
        )
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001, 0.0, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 133.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();

    // Use a simulation bound so that atoms that escape the capture region are deleted from the simulation - glass box.
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, shift_z * 0.001),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.05,  0.02, 0.02),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // The simulation bound also now includes a narrow pipe connecting to the 3D area MOT.
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.06, 0.0, shift_z * 0.001),
        })
        .with(Cylinder {
            radius: 0.00075,
            length: 0.02,
            direction: dir,
            perp_x: perp_x,
            perp_y: perp_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Adding a tube leading to the 3D MOT
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.5*(distance * 0.001 - 0.07) + 0.07, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.5*(distance * 0.001 - 0.07), 0.02, 0.02),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    // Adding a space of the 3D MOT
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(distance * 0.001, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(0.05, 0.05, 0.05),
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

    // sim.world.insert(EmissionForceOption::On(EmissionForceConfiguration {
    //     explicit_threshold: 5,
    // }));
    // sim.world.insert(ScatteringFluctuationsOption::On);

    // Run the simulation for a number of steps.
    for _i in 0..(time * 1000) {
        sim.step();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}