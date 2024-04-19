//
//  Particle.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 4/12/24.
//

import Foundation
import Observation

@Observable class Particle1D{
    var position: Double
    var velocity: Double
    var force: [Double]
    var imagePositions: [Double]
    let boxSize: Double
    
    init(position: Double, velocity: Double, boxSize: Double) {
        self.position = position
        self.velocity = velocity
        self.force = Array(repeating: 0, count: 2)
        self.boxSize = boxSize
        self.imagePositions = []
        self.updateImages()
    }
    
    func updateImages() {   // 1D
        self.imagePositions = []
        for i in stride(from: -1, through: 1, by: 1) {
            self.imagePositions.append(self.position + Double(i)*boxSize)
        }
    }
    
    func changePosition(newPosition: Double){
        self.position = newPosition
        self.updateImages()
    }
}
