extern crate atomecs as lib;
extern crate nalgebra;
use lib::atom::{Atom, Position, Velocity};
use lib::atom_sources::AtomSourcePlugin;
use lib::integrator::Timestep;
use lib::laser::LaserPlugin;
use lib::laser::gaussian::GaussianBeam;
// use lib::laser::gaussian::CircularMask;
use lib::laser_cooling::{CoolingLight, LaserCoolingPlugin};
use lib::magnetic::quadrupole::QuadrupoleField3D;
use lib::output::file::FileOutputPlugin;
use lib::output::file::Text;
use lib::shapes::{Cuboid, Cylinder};
// use lib::gravity::ApplyGravitationalForceSystem;
// use lib::gravity::ApplyGravityOption;
use lib::sim_region::{SimulationVolume, VolumeType};
use lib::simulation::SimulationBuilder;
use lib::species::{Potassium39, Potassium39_767D2};
use lib::atom_sources::oven::{OvenAperture, OvenBuilder};
use lib::atom_sources::mass::{MassDistribution, MassRatio};
use lib::atom_sources::VelocityCap;
use lib::atom_sources::emit::AtomNumberToEmit;
use lib::destructor::ToBeDestroyed;
use nalgebra::Vector3;
use specs::prelude::*;
use std::time::Instant;
extern crate specs;
#[allow(unused_imports)]
use specs::{Builder, Entity, RunNow, World};
#[allow(unused_imports)]

const BEAM_NUMBER : usize = 10;

pub fn calculate_rayleigh_range(wavelength: &f64, e_radius: &f64) -> f64 {
    2.0 * 3.14 * e_radius.powf(2.0) / wavelength
}

fn main() {
    let now = Instant::now();

    let mut sim_builder = SimulationBuilder::default();
    sim_builder.add_plugin(LaserPlugin::<{BEAM_NUMBER}>);
    sim_builder.add_plugin(LaserCoolingPlugin::<Potassium39_767D2, {BEAM_NUMBER}>::default());
    sim_builder.add_plugin(AtomSourcePlugin::<Potassium39>::default());
    sim_builder.add_plugin(FileOutputPlugin::<Position, Text, Atom>::new("pos.txt".to_string(), 10));
    sim_builder.add_plugin(FileOutputPlugin::<Velocity, Text, Atom>::new("vel.txt".to_string(), 10));
    let mut sim = sim_builder.build();

    // Shift in z axis between the 2D and 3D MOT and the x axis distance between them 'It equals 347 mm and includes a 55.7 mm-long differential pumping section, which is inserted into the 2D-MOT chamber'
    let shift = 0.0;
    let thick = 1.0;

    //I changed 1.4142*small_radius to b small_radius because I dont know what this 1.4142 is

    // Magnetic gradient
    let magnetic_grad_2d = 17.0;

    // Push beam along x axis
    let detuning_push = -5.2 * 6.2;
    let power_push = 6.0;
    let radius_push = 0.75;
    let wavelenght = 767.5e-9;
    let radius_pipe = 0.5;

    // // Create cooling lasers - 2DMOT.
    let detuning_2d = -5.8 * 6.2;
    let power_2d = 80.0;
    let radius_2d = 14.1;
    let ellipticity = 0.943;

    let number_to_emit = 10000000.0;
    let velocity_cap = 250.0;
    let temperature = 316.15;

    let time = 6;
    let small_radius = ((-1.0 * (ellipticity * ellipticity * radius_2d * 0.001 * radius_2d * 0.001 - radius_2d * 0.001 * radius_2d * 0.001)) as f64).sqrt();

    let number_to_emit_b = (number_to_emit*(small_radius*2.0*small_radius*2.0/(2.0*small_radius*2.0*small_radius*2.0 + 4.0*small_radius*2.0*radius_2d * 0.001*2.0))) as i32;
    let number_to_emit_s = (number_to_emit*(small_radius*2.0*radius_2d * 0.001*2.0/(2.0*small_radius*2.0*small_radius*2.0 + 4.0*small_radius*2.0*radius_2d * 0.001*2.0))) as i32;

    let dir = Vector3::new(1.0, 0.0, 0.0).normalize();
    let perp_x = dir.normalize().cross(&dir);
    let perp_y = dir.normalize().cross(&perp_x);

    // sim.world.insert(ApplyGravityOption);

    //Create magnetic field 2D MOT.
    sim.world
        .create_entity()
        .with(QuadrupoleField3D::gauss_per_cm(magnetic_grad_2d, Vector3::x()))
        .with(Position {pos: Vector3::new(0.0, 0.0, shift * 0.001)})
        .build();

    // Push beam
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift * 0.001),
            e_radius: radius_push * 0.001 / 2.0_f64.sqrt(),
            power: power_push * 0.001,
            direction: Vector3::new(1.0,0.0,0.0).normalize(),
            rayleigh_range: calculate_rayleigh_range(&wavelenght, &(radius_push * 0.001)),
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_push,
            -1,
        ))
        // .with(CircularMask{
        //     radius: radius_pipe * 0.001,
        // })
        .build();

    // 4 cooling lasers of 2D MOT cooling
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001,
            direction: Vector3::new(0.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001,
            direction: Vector3::new(0.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001,
            direction: Vector3::new(0.0, 0.0, -1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_2d,
            1,
        ))
        .build();
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(0.0, 0.0, shift * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * 0.001,
            direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_2d,
            1,
        ))
        .build();

    // 4 cooling lasers of 2D MOT repump
    // sim.world
    //     .create_entity()
    //     .with(GaussianBeam {
    //         intersection: Vector3::new(0.0, 0.0, shift * 0.001),
    //         e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
    //         power: 50.0 * 0.001,
    //         direction: Vector3::new(0.0, 1.0, 0.0).normalize(),
    //         rayleigh_range: f64::INFINITY,
    //         ellipticity: ellipticity,
    //     })
    //     .with(CoolingLight::for_transition::<Potassium39_767D2>(
    //         -3.9 * 6.2,
    //         1,
    //     ))
    //     .build();
    // sim.world
    //     .create_entity()
    //     .with(GaussianBeam {
    //         intersection: Vector3::new(0.0, 0.0, shift * 0.001),
    //         e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
    //         power: 50.0 * 0.001,
    //         direction: Vector3::new(0.0, -1.0, 0.0).normalize(),
    //         rayleigh_range: f64::INFINITY,
    //         ellipticity: ellipticity,
    //     })
    //     .with(CoolingLight::for_transition::<Potassium39_767D2>(
    //         -3.9 * 6.2,
    //         1,
    //     ))
    //     .build();
    // sim.world
    //     .create_entity()
    //     .with(GaussianBeam {
    //         intersection: Vector3::new(0.0, 0.0, shift * 0.001),
    //         e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
    //         power: 50.0 * 0.001,
    //         direction: Vector3::new(0.0, 0.0, -1.0).normalize(),
    //         rayleigh_range: f64::INFINITY,
    //         ellipticity: ellipticity,
    //     })
    //     .with(CoolingLight::for_transition::<Potassium39_767D2>(
    //         -3.9 * 6.2,
    //         1,
    //     ))
    //     .build();
    // sim.world
    //     .create_entity()
    //     .with(GaussianBeam {
    //         intersection: Vector3::new(0.0, 0.0, shift * 0.001),
    //         e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
    //         power: 50.0 * 0.001,
    //         direction: Vector3::new(0.0, 0.0, 1.0).normalize(),
    //         rayleigh_range: f64::INFINITY,
    //         ellipticity: ellipticity,
    //     })
    //     .with(CoolingLight::for_transition::<Potassium39_767D2>(
    //         -3.9 * 6.2,
    //         1,
    //     ))
    //     .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - upper big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(0.0, 0.0, -1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001*2.0,small_radius*2.0, thick * 0.001]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, 0.0, small_radius),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - lower big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(0.0, 0.0, 1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001*2.0,small_radius*2.0, thick * 0.001]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, 0.0, -small_radius),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - back big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(0.0, -1.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001*2.0, thick * 0.001,small_radius*2.0]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, small_radius, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - front big plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(0.0, 1.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001*2.0, thick * 0.001,small_radius*2.0]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(0.0, -small_radius, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_b,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - left small plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(1.0, 0.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [ thick * 0.001,small_radius*2.0,small_radius*2.0]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(-radius_2d * 0.001, 0.0, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();
    
    // Atoms are taken from background gas which is inserted by 6 'ovens' - right small plate
    // Add atoms
    sim.world
        .create_entity()
        .with(
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(-1.0, 0.0, 0.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [thick * 0.001,small_radius*2.0,small_radius*2.0]})
                .build(),
        )
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001, 0.0, 0.0),
        })
        .with(MassDistribution::new(vec![MassRatio {
            mass: 39.0,
            ratio: 1.0,
        }]))
        // .with(EmitNumberPerFrame { number: 28 })
        .with(AtomNumberToEmit {
            number: number_to_emit_s,
        })
        .with(ToBeDestroyed)
        .build();

    // 2D MOT chamber
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, shift * 0.001),
        })
        .with(Cuboid {
            half_width: Vector3::new(radius_2d * 0.001 + thick * 0.001*0.5, small_radius + thick * 0.001*0.5, small_radius + thick * 0.001*0.5),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //Small hole
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001 + thick * 0.001*0.5 + 0.01*0.5, 0.0, shift * 0.001),
        })
        .with(Cylinder {
            radius: radius_pipe*0.001,
            length: 0.01,
            direction: dir,
            perp_x: perp_x,
            perp_y: perp_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //Some propagation
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001 + thick * 0.001*0.5 + 0.01 + 0.005, 0.0, shift * 0.001),
        })
        .with(Cylinder {
            radius: 0.0047,
            length: 0.01,
            direction: dir,
            perp_x: perp_x,
            perp_y: perp_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //My tube extracting beam at 34 mrad cone
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001 + thick * 0.001*0.5 + 0.01 + 0.01 + 0.0005, 0.0, shift * 0.001),
        })
        .with(Cylinder {
            radius: 0.00117685,
            length: 0.001,
            direction: dir,
            perp_x: perp_x,
            perp_y: perp_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    //My tube allowing to propagate
    sim.world
        .create_entity()
        .with(Position {
            pos: Vector3::new(radius_2d * 0.001 + thick * 0.001*0.5 + 0.01 + 0.01 + 0.001 + 0.15, 0.0, shift * 0.001),
        })
        .with(Cylinder {
            radius: 0.1,
            length: 0.3,
            direction: dir,
            perp_x: perp_x,
            perp_y: perp_y
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

    sim.world.insert(VelocityCap { value: velocity_cap });

    // Define timestep
    sim.world.insert(Timestep { delta: 1.0e-6 });

    // sim.world.insert(ApplyGravitationalForceSystem);

    // Run the simulation for a number of steps.
    for _i in 0..(time * 1000) {
        sim.step();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}