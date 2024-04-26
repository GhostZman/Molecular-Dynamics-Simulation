//
//  MD Simulation.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

import Foundation
import Observation
import SwiftUI

@Observable class MDSimulation: ObservableObject{
    var boxLength: Int = 1
    var nStep: Int = Int(1E4)
    var stepSize: Double = 0.004
    var halfTimeStep: Double = 0.002
    var PE: Double = 0.0
    var KE: Double = 0.0
    var temperature: Double = 0.0
    var initialTemperature: Double = 0.0
    let nDim: Int = 3
    let nPrint: Int = 100
    var conditionsSet: Bool = false
    var cubrtUnitCells: Int = 2
    
    var atoms: [Particle3D] = []

    
    //let objectWillChange = PassthroughSubject<Void, Never>()
    
    func randVelocityVector() -> [Double]{ // Inital velocities
        var atomVelocity = (Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0)+Double.random(in: 0.0 ... 1.0))/12.0 - 0.5
        let phiDirection = Double.random(in: 0.0 ... Double.pi)
        let thetaDirection = Double.random(in:0.0 ... Double.pi)
        atomVelocity = atomVelocity*sqrt(initialTemperature)   // Scale v with temperature
        let velocityVector = [atomVelocity*cos(phiDirection)*sin(thetaDirection), atomVelocity*sin(phiDirection)*sin(thetaDirection), atomVelocity*cos(thetaDirection)]
        return velocityVector
    }
    
    func setInitialTemperature(temp: Double){
        self.initialTemperature = temp
    }
    
    func makeParticles(option: String) {
        
        let FCCPoints: [[Double]] = [[0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]]
        
        switch option {
        case "FCC Lattice":
            self.boxLength = 2
            for pos in FCCPoints {
                atoms.append(Particle3D(position: pos, velocity: randVelocityVector(),  boxSize: Double(boxLength)))
            }
        case "FCC Lattice 2":
            self.boxLength = 2 * cubrtUnitCells
            for i in stride(from: 0, to: cubrtUnitCells, by: 1) {
                for j in stride(from: 0, to: cubrtUnitCells, by: 1) {
                    for k in stride(from: 0, to: cubrtUnitCells, by: 1) {
                        for pos in FCCPoints {
                            atoms.append(Particle3D(position: [pos[0] + 2*Double(i), pos[1] + 2*Double(j), pos[2] + 2*Double(k)], velocity: randVelocityVector(),  boxSize: Double(boxLength)))
                        }
                    }
                }
            }
            
        default:
            var atomIndex: Int
            
            self.boxLength = 4
            
            
            atomIndex = -1
            for ix in stride(from: 0, through: boxLength-1, by: 1){  //Set up lattice of side boxLength
                atomIndex = atomIndex+1
                let atomPosition = [Double(ix), Double(ix), Double(ix)]
                atoms.append(Particle3D(position: atomPosition, velocity: randVelocityVector(),  boxSize: Double(boxLength)))
                
                print("init atomVelocity = \(atoms[atomIndex].velocity)")
                
            }
            print("numAtoms = \(atoms.count) boxLength = \(boxLength)")
        }
        self.conditionsSet = true
    }
    
    func runSimulation() {
        if conditionsSet {
            
            var timeIndex: Int
            
            timeIndex = 0
            KE = 0.0    // inital KE & PE
            PE = 0.0
            PE = updateForces()
            for atomIndex in stride(from: 0, through: atoms.count - 1, by: 1) {
                KE = KE + (pow(atoms[atomIndex].velocityMagnitude(), 2))/2
            }
            //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
            print("step: \(timeIndex)\nposition    velocity         force ")
            atoms.forEach { particle in
                print(particle.position, particle.velocity, particle.force)
            }
            
            for timeIndex in stride(from: 1, to: nStep, by: 1){     // Main loop
                //usleep(5000)
                //objectWillChange.send()
                for atomIndex in stride(from: 0, through: atoms.count - 1, by: 1) {   // Velocity Verlet
                    for dim in stride(from: 0, through: 2, by: 1){
                        var newPosition = atoms[atomIndex].position[dim] + stepSize*(atoms[atomIndex].velocity[dim] + halfTimeStep*atoms[atomIndex].force[dim][1])      //PBC
                        if (newPosition <= 0.0){
                            newPosition = newPosition + Double(boxLength)
                        }
                        if (newPosition >= Double(boxLength)){
                            newPosition = newPosition - Double(boxLength)
                        }
                        atoms[atomIndex].changePosition(newPosition: newPosition, dimIndex: dim)
                    }
                }
                PE = updateForces()
                KE = 0
                for atomIndex in stride(from: 0, through: atoms.count - 1, by: 1){
                    for dim in stride(from: 0, through: 1, by: 1){
                        atoms[atomIndex].velocity[dim] = atoms[atomIndex].velocity[dim] + halfTimeStep*(atoms[atomIndex].force[dim][1] + atoms[atomIndex].force[dim][0])
                    }
                    KE = KE + (atoms[atomIndex].velocityMagnitude()*atoms[atomIndex].velocityMagnitude())/2
                }
                temperature = 2.0 * KE / (3.0 * Double(atoms.count))
                print("step: \(timeIndex)\nposition")
                atoms.forEach { particle in
                    print(particle.position)
                }
                if (timeIndex % nPrint == 0){
                    //print("\(timeIndex) PE = \(PE) KE = \(KE) PE+KE = \(PE+KE)")
                    /*print("step: \(timeIndex)\nposition    velocity         force ")
                     atoms.forEach { particle in
                     print(particle.position, particle.velocity, particle.force)
                     }
                     */
                }
            }
        } else {
            print("Inital conditions not set")
        }
    }
    
    func updateForces() -> Double {   //run this AFTER updating positions BEFORE updating velocities
        var forceOfiOnj: Double
        var separation2: Double
        var inverseSeparation2: Double = 0.0
        var displacementX: Double
        var displacementY: Double
        var displacementZ: Double
        var phiAngle: Double
        var thetaAngle: Double
        let cutoffSeparation2: Double = 9.0
        
        var newPE: Double = 0.0     // Initialize
        for i in stride(from: 0, through: atoms.count - 1, by: 1) {
            for dim in stride(from: 0, through: 2, by: 1) {
                atoms[i].force[dim][0] = atoms[i].force[dim][1]
                atoms[i].force[dim][1] = 0
            }
        }
        
        for i in stride(from: 0, through: atoms.count - 2, by: 1){
            for j in stride(from: i+1, through: atoms.count - 1, by: 1) {
                for atomImage in atoms[j].imagePositions{
                    displacementX = atoms[i].position[0] - atomImage[0]
                    displacementY = atoms[i].position[1] - atomImage[1]
                    displacementZ = atoms[i].position[2] - atomImage[2]
                    
                    thetaAngle = acos(displacementZ/sqrt(pow(displacementX, 2)+pow(displacementY, 2)+pow(displacementZ, 2)))
                    phiAngle = atan2(displacementY, displacementX)
                    separation2 = pow(displacementX, 2) + pow(displacementY, 2) + pow(displacementZ, 2)
                    if (separation2 < cutoffSeparation2) {   // Cut off
                        if (abs(separation2) < 1.0E-7) {
                            separation2 = 1.0E-7
                        }
                        inverseSeparation2 = 1.0/separation2
                        forceOfiOnj = 48*(pow(inverseSeparation2, 3.0) - 0.5)*pow(inverseSeparation2, 3.0)
                        forceOfiOnj = forceOfiOnj*inverseSeparation2*sqrt(separation2)
                        
                        atoms[i].force[0][1] = atoms[i].force[0][1] + forceOfiOnj*sin(thetaAngle)*cos(phiAngle)
                        atoms[j].force[0][1] = atoms[j].force[0][1] - forceOfiOnj*sin(thetaAngle)*cos(phiAngle)
                        
                        atoms[i].force[1][1] = atoms[i].force[1][1] + forceOfiOnj*sin(thetaAngle)*sin(phiAngle)
                        atoms[j].force[1][1] = atoms[j].force[1][1] - forceOfiOnj*sin(thetaAngle)*sin(phiAngle)
                        
                        atoms[i].force[2][1] = atoms[i].force[2][1] + forceOfiOnj*cos(thetaAngle)
                        atoms[j].force[2][1] = atoms[j].force[2][1] - forceOfiOnj*cos(thetaAngle)
                        
                        newPE = newPE + 4 * pow(inverseSeparation2, 3) * (pow(inverseSeparation2, 3) - 1)
                    }
                }
            }
        }
        return newPE
    }
    func sign(a: Double, b: Double) -> Double{
        if (b >= 0) {
            return abs(a)
        }
        else {
            return -abs(a)
        }
    }
    
}
