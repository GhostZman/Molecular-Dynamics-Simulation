//
//  DrawingView.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 4/19/24.
//
import SwiftUI

struct DrawingView: View {
    var atoms: [Particle2D]
    
    var body: some View {
        GeometryReader { geometry in
            ZStack {
                ForEach(atoms.indices, id: \.self) { index in
                    Circle()
                        .fill(Color.blue)
                        .frame(width: 10, height: 10)
                        .position(x: CGFloat(atoms[index].position[0]*25), y: CGFloat(atoms[index].position[1]*25))
                }
            }
            .frame(width: geometry.size.width, height: geometry.size.height)
        }
    }
}
