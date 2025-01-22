//
//  MoleculeView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SceneKit

struct MoleculeView: UIViewRepresentable {
    func makeUIView(context: Context) -> SCNView {
        let sceneView = SCNView()
        
        sceneView.scene = createMoleculeScene()
        sceneView.allowsCameraControl = true
        sceneView.backgroundColor = .black
        
        return sceneView
        
    }
    
    func updateUIView(_ uiView: SCNView, context: Context) {}
    
    private func createMoleculeScene() -> SCNScene {
        let scene = SCNScene()
        
        let oxygen = SCNSphere(radius: 0.3)
        oxygen.firstMaterial?.diffuse.contents = UIColor.red
        let oxygenNode = SCNNode(geometry: oxygen)
        scene.rootNode.addChildNode(oxygenNode)
        
        let hydrogen1 = SCNSphere(radius: 0.2)
        hydrogen1.firstMaterial?.diffuse.contents = UIColor.white
        
        let hydrogenNode1 = SCNNode(geometry: hydrogen1)
        hydrogenNode1.position = SCNVector3(-0.8, 0, 0)
        scene.rootNode.addChildNode(hydrogenNode1)
        
        let hydrogen2 = SCNSphere(radius: 0.2)
        hydrogen2.firstMaterial?.diffuse.contents = UIColor.white
        
        let hydrogenNode2 = SCNNode(geometry: hydrogen2)
        hydrogenNode2.position = SCNVector3(0.8, 0, 0)
        scene.rootNode.addChildNode(hydrogenNode2)
        
        let bond1 = SCNCylinder(radius: 0.05, height: 0.8)
        bond1.firstMaterial?.diffuse.contents = UIColor.gray
        
        let bondNode1 = SCNNode(geometry: bond1)
        bondNode1.position = SCNVector3(-0.4, 0, 0)
        bondNode1.eulerAngles = SCNVector3(0, 0, (Double.pi / 2))
        scene.rootNode.addChildNode(bondNode1)
        
        let bond2 = SCNCylinder(radius: 0.05, height: 0.8)
        bond2.firstMaterial?.diffuse.contents = UIColor.gray
        
        let bondNode2 = SCNNode(geometry: bond2)
        bondNode2.position = SCNVector3(0.4, 0, 0)
        bondNode2.eulerAngles = SCNVector3(0, 0, (Double.pi / 2))
        scene.rootNode.addChildNode(bondNode2)
        
        return scene
        
    }
}

#Preview {
    MoleculeView()
        .ignoresSafeArea(edges: .all)
    
}
