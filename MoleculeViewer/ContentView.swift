//
//  ContentView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SwiftData
import UIKit

struct ContentView: View {
    let scenes: [String] = ["Structures", "Electronegativity", "Shapes", "Polarity", "Simulation"]
    @State var selectedScene: String = "Structures"
    
    var body: some View {
        ZStack {
            Rectangle()
                .fill(.red)
                .frame(width: 200, height: 200)
            
            MoleculeSimulationView()
            
            VStack {
                Spacer()
                
                ZStack {
                    Capsule()
                        .fill(.white)
                        .opacity(0.5)
                        .blur(radius: 30)
                        .frame(width: 250, height: 100)
                    
                    VStack {
//                    Text("Selection: \(selectedScene)")
//                        .foregroundStyle(.black)
                    
                    Text("N\u{2082}")
                        .font(.system(size: 32, weight: .bold, design: .default))
                        .foregroundStyle(.black)
                        
                    Text("Nitrogen Gas")
                        .font(.system(size: 24, weight: .medium, design: .default))
                        .foregroundStyle(.black)
                    
//                        Text("C\u{2088}H\u{2081}\u{2080}N\u{2084}O\u{2082}")
//                            .font(.system(size: 32, weight: .bold, design: .default))
//                            .foregroundStyle(.black)
//                        
//                        Text("Caffeine")
//                            .font(.system(size: 24, weight: .medium, design: .default))
//                            .foregroundStyle(.black)
                    }
                }
            }
        }
        .statusBar(hidden: true)
        
    }
}

#Preview {
    ContentView()
    
}
