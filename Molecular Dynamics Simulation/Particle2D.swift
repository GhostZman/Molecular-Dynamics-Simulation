//
//  Particle2D.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 4/19/24.
//

import Foundation
import Observation

@Observable class Particle2D{
    var position: [Double]
    var velocity: [Double]
    var force: [[Double]]
    var imagePositions: [[Double]]
    let boxSize: Double
    
    init(position: [Double], velocity: [Double], boxSize: Double) {
        if position.count == 2 {
            self.position = position
        } else {
            self.position = Array(repeating: 0, count: 2)
            print("Position has wrong number of dimensions, setting to 0")
        }
        if velocity.count == 2 {
            self.velocity = velocity
        } else {
            self.velocity = Array(repeating: 0, count: 2)
            print("Velocity has wrong number of dimensions, setting to 0")
        }
        self.force = Array(repeating: Array(repeating: 0, count: 2), count: 2)
        self.boxSize = boxSize
        self.imagePositions = []
        self.updateImages()
    }
    
    func updateImages() {
        self.imagePositions = []
        for i in stride(from: -1, through: 1, by: 1) {
            for j in stride(from: -1, through: 1, by: 1){
                self.imagePositions.append([self.position[0] + Double(i)*boxSize, self.position[1] + Double(j)*boxSize])
            }
        }
    }
    
    func changePosition(newPosition: [Double]) {
        if newPosition.count == 2 {
            self.position = newPosition
        } else {
            print("new position has wrong number of dimensions, it will not be changed")
        }
        self.updateImages()
    }
    
    func changePosition(newPosition: Double, dimIndex: Int){
        self.position[dimIndex] = newPosition
        self.updateImages()
    }
    
    func velocityMagnitude() -> Double {
        return sqrt(pow(velocity[0], 2) + pow(velocity[1], 2))
    }
}
