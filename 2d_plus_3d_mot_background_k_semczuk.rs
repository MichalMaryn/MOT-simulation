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
use lib::laser_cooling::{CoolingLight, LaserCoolingPlugin};
use lib::magnetic::quadrupole::QuadrupoleField3D;
use lib::output::file::FileOutputPlugin;
use lib::output::file::Text;
use lib::shapes::{Cuboid, Cylinder};
use lib::gravity::ApplyGravitationalForceSystem;
use lib::gravity::ApplyGravityOption;
use lib::sim_region::{SimulationVolume, VolumeType};
use lib::simulation::SimulationBuilder;
use lib::species::{Potassium39, Potassium39_767D2};
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
    sim_builder.add_plugin(LaserCoolingPlugin::<Potassium39_767D2, {BEAM_NUMBER}>::default());
    sim_builder.add_plugin(AtomSourcePlugin::<Potassium39>::default());
    sim_builder.add_plugin(FileOutputPlugin::<Position, Text, Atom>::new("pos.txt".to_string(), 100));
    sim_builder.add_plugin(FileOutputPlugin::<Velocity, Text, Atom>::new("vel.txt".to_string(), 100));
    // sim_builder.add_plugin(FileOutputPlugin::<MagneticFieldSampler, Text, Atom>::new("mag.txt".to_string(), 1));
    let mut sim = sim_builder.build();

    let shift_x = 0.0;
    let shift_z = 1.5;
    let thick = 1.0;

    sim.world.insert(ApplyGravityOption);

    // Magnetic gradient
    let quadrupole_gradient_2d = 3.938;
    let quadrupole_gradient_3d = 8.5; //23.708;

    // Push beam along x axis
    let detuning_push= -32.988; //-30.042
    let power_push = 0.51383;
    let radius_push = 1.95;
    let wavelenght= 767.0e-9;

    // Best observed small 3d 22/05: power2d 195.71   power3d 98.106   detuning -40.487         gradient3d 16.428        shiftz 50.383 
    // Best observed 02/06: 
    // radius_2d    detuning    ellipticity    ratio_power_3d_2d    quadrupole_gradient_2d    quadrupole_gradient_3d
    // _________    ________    ___________    _________________    ______________________    ______________________
    //  17.514      -31.714       91.571            69.014                  5.5191                    11.864

    // 03/06
    // Best estimated feasible point (according to models):
    //     detuning    power_push    radius_push    radius_2d    ellipticity    ratio_power_3d_2d    quadrupole_gradient_2d    quadrupole_gradient_3d
    //     ________    __________    ___________    _________    ___________    _________________    ______________________    ______________________

    //     -32.988       51.383        195.45        18.398        94.073            51.282                  3.938                     23.708 
    // Estimated objective function value = -2.546e-05

    // // Create cooling lasers - 2DMOT
    let detuning_2d = -32.988;
    let power_2d = 300.0;
    let radius_2d = 18.398;
    let ellipticity = 0.94073;

    let detuning_3d = -32.988;
    // let power_3d = 98.106;
    let radius_3d = 7.5;

    let ratio_power_3d_2d = 0.51282;

    let number_to_emit = 10000000.0;
    let velocity_cap = 70.0;
    let temperature = 316.0;

    let distance = 350.0;

    let time = 60;
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
        .with(Position {pos: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001)})
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
            intersection: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001),
            e_radius: radius_push * 0.001 / 2.0_f64.sqrt(),
            power: power_push * 0.001,
            direction: Vector3::new(1.0,0.0,0.0).normalize(),
            rayleigh_range: calculate_rayleigh_range(&wavelenght, &radius_push),
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
            detuning_push,
            -1,
        ))
        .build();

    // Add 4 2D MOT lasers
    sim.world
        .create_entity()
        .with(GaussianBeam {
            intersection: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 1.0, 1.0).normalize(),
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
            intersection: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, -1.0, -1.0).normalize(),
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
            intersection: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, 1.0, -1.0).normalize(),
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
            intersection: Vector3::new(shift_x * 0.001, 0.0, shift_z * 0.001),
            e_radius: radius_2d * 0.001 / 2.0_f64.sqrt(),
            power: power_2d * ratio_power_3d_2d * 0.001 / 2.0,
            direction: Vector3::new(0.0, -1.0, 1.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: ellipticity,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
            power: power_2d * (1.0 - ratio_power_3d_2d) * 0.001 / 3.0,
            direction: Vector3::new(1.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
            direction: Vector3::new(-1.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
            direction: Vector3::new(-1.0, 1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
            direction: Vector3::new(1.0, -1.0, 0.0).normalize(),
            rayleigh_range: f64::INFINITY,
            ellipticity: 0.0,
        })
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
        .with(CoolingLight::for_transition::<Potassium39_767D2>(
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
            OvenBuilder::<Potassium39>::new(
                temperature,
                Vector3::new(0.0, 0.0, -1.0).normalize())
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, small_radius*2.0, thick * 0.001]})
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
            number: number_to_emit_s,
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
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, small_radius*2.0, thick * 0.001]})
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
            number: number_to_emit_s,
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
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, thick * 0.001, small_radius*2.0]})
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
            number: number_to_emit_s,
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
                .with_aperture(OvenAperture::Cubic { size: [radius_2d * 0.001 * 2.0, thick * 0.001, small_radius*2.0]})
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
            number: number_to_emit_s,
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
                .with_aperture(OvenAperture::Cubic { size: [thick * 0.001, small_radius * 2.0, small_radius * 2.0]})
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
            number: number_to_emit_b,
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
                .with_aperture(OvenAperture::Cubic { size: [thick * 0.001, small_radius*2.0, small_radius*2.0]})
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

    // Run the simulation for a number of steps.
    for _i in 0..(time * 1000) {
        sim.step();
    }

    // sim.world.insert(ApplyCollisionsOption);
    // sim.world.insert(CollisionParameters {
    //     macroparticle: 4e2,
    //     box_number: 200, //Any number large enough to cover entire cloud with collision boxes. Overestimating box number will not affect performance.
    //     box_width: 20e-6, //Too few particles per box will both underestimate collision rate and cause large statistical fluctuations.
    //     //Boxes must also be smaller than typical length scale of density variations within the cloud, since the collisions model treats gas within a box as homogeneous.
    //     sigma: 3.5e-16, //Approximate collisional cross section of Rb87
    //     collision_limit: 10_000_000.0, //Maximum number of collisions that can be calculated in one frame.
    //                                    //This avoids absurdly high collision numbers if many atoms are initialised with the same position, for example.
    // });
    // sim.world.insert(CollisionsTracker {
    //     num_collisions: Vec::new(),
    //     num_atoms: Vec::new(),
    //     num_particles: Vec::new(),
    // });

    // sim.world.insert(EmissionForceOption::On(EmissionForceConfiguration {
    //     explicit_threshold: 5,
    // }));
    // sim.world.insert(ScatteringFluctuationsOption::On);

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}