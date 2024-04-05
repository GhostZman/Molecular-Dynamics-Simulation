//
//  ContentView.swift
//  Molecular Dynamics Simulation
//
//  Created by Phys440Zachary on 3/29/24.
//

import SwiftUI

struct ContentView: View {
    var body: some View {
        VStack {
            Canvas { context, size in
                context.stroke(Path(roundedRect: CGRect(origin: .zero, size: size), cornerSize: CGSize(width: 20, height: 50)), with: .color(.red))
                context.stroke(Path(roundedRect: CGRect(origin: .zero, size: size), cornerSize: CGSize(width: 50, height: 20)), with: .color(.green))
            }
            .frame(width: 300, height: 200)
        }
        .background(.black)
        .ignoresSafeArea()
        .padding()
    }
}
