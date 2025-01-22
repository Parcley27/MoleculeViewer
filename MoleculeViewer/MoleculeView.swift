//
//  MoleculeView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SceneKit

struct MoleculeView: UIViewRepresentable {
    let molFile: String
    
    func makeUIView(context: Context) -> SCNView {
        let sceneView = SCNView()
        sceneView.scene = createMoleculeScene(from: molFile)
        sceneView.allowsCameraControl = true
        sceneView.backgroundColor = .black
        return sceneView
    }

    func updateUIView(_ uiView: SCNView, context: Context) {}

    private func createMoleculeScene(from molFile: String) -> SCNScene {
        let scene = SCNScene()
        let lines = molFile.components(separatedBy: .newlines)
        
        // Read atom and bond counts from the 4th line
        let countsLine = lines[3]
        let atomCountRange = countsLine.index(countsLine.startIndex, offsetBy: 0)..<countsLine.index(countsLine.startIndex, offsetBy: 3)
        let bondCountRange = countsLine.index(countsLine.startIndex, offsetBy: 3)..<countsLine.index(countsLine.startIndex, offsetBy: 6)
        
        
        guard let atomCount = Int(countsLine[atomCountRange].trimmingCharacters(in: .whitespaces)),
              let bondCount = Int(countsLine[bondCountRange].trimmingCharacters(in: .whitespaces)) else {
            return scene
            
        }
        
        // Parse atoms
        var atomPositions = [SCNVector3]()
        for i in 0..<atomCount {
            let atomLine = lines[4 + i]
            let components = atomLine.split(separator: " ", omittingEmptySubsequences: true)
            if components.count >= 4 {
                let x = Float(components[0]) ?? 0
                let y = Float(components[1]) ?? 0
                let z = Float(components[2]) ?? 0
                atomPositions.append(SCNVector3(x, y, z))
                
                // Create and add atom node
                let atomNode = createAtom(radius: 0.2, color: .cyan)
                atomNode.position = SCNVector3(x, y, z)
                scene.rootNode.addChildNode(atomNode)
            }
        }
        
        // Parse bonds
        for i in 0..<bondCount {
            let bondLine = lines[4 + atomCount + i]
            let components = bondLine.split(separator: " ", omittingEmptySubsequences: true)
            if components.count >= 3 {
                let startIndex = Int(components[0]) ?? 0
                let endIndex = Int(components[1]) ?? 0
                
                if startIndex > 0 && endIndex > 0 {
                    let start = atomPositions[startIndex - 1]
                    let end = atomPositions[endIndex - 1]
                    
                    // Create and add bond node
                    let bondNode = createBond(from: start, to: end, radius: 0.05, color: .gray)
                    scene.rootNode.addChildNode(bondNode)
                }
            }
        }
        
        return scene
    }
    
    private func createAtom(radius: CGFloat, color: UIColor) -> SCNNode {
        let sphere = SCNSphere(radius: radius)
        sphere.firstMaterial?.diffuse.contents = color
        return SCNNode(geometry: sphere)
    }
    
    private func createBond(from start: SCNVector3, to end: SCNVector3, radius: CGFloat, color: UIColor) -> SCNNode {
        let vector = SCNVector3(end.x - start.x, end.y - start.y, end.z - start.z)
        let height = sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z)
        
        let cylinder = SCNCylinder(radius: radius, height: CGFloat(height))
        cylinder.firstMaterial?.diffuse.contents = color
        
        let bondNode = SCNNode(geometry: cylinder)
        bondNode.position = SCNVector3((start.x + end.x) / 2, (start.y + end.y) / 2, (start.z + end.z) / 2)
        bondNode.look(at: end)
        bondNode.eulerAngles.x += .pi / 2
        
        return bondNode
        
    }
}

#Preview {
    let testData = """
5793
  -OEChem-01222501053D

 24 24  0     1  0  0  0  0  0999 V2000
   -0.6679    1.1587    0.2570 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8870   -2.4483   -0.3388 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8623   -2.0693    0.4696 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8609    0.5414   -0.4619 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1222    2.6552    0.2574 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3742    0.9717   -0.1865 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727   -1.2470    0.2300 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.0856   -1.0709   -0.1940 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.2211   -0.0621   -0.2375 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.6082    0.3151    0.1839 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.6388    1.4132   -0.2534 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.6550   -0.1577    0.2740 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4248   -1.3522    1.3206 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2066   -1.2487   -1.2697 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2548   -0.0098   -1.3343 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7952    0.3598    1.2636 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5967    1.5141   -1.3440 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6916   -0.1535    1.3685 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1564   -1.0581   -0.0922 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8514   -2.3615   -1.3066 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4973   -2.9356    0.2200 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7165    0.4989   -1.4227 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4876    2.5033    1.1448 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9192    1.7652    0.1440 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  1  0  0  0  0
  1 11  1  0  0  0  0
  2  7  1  0  0  0  0
  2 20  1  0  0  0  0
  3  8  1  0  0  0  0
  3 21  1  0  0  0  0
  4 10  1  0  0  0  0
  4 22  1  0  0  0  0
  5 11  1  0  0  0  0
  5 23  1  0  0  0  0
  6 12  1  0  0  0  0
  6 24  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 13  1  0  0  0  0
  8 10  1  0  0  0  0
  8 14  1  0  0  0  0
  9 12  1  0  0  0  0
  9 15  1  0  0  0  0
 10 11  1  0  0  0  0
 10 16  1  0  0  0  0
 11 17  1  0  0  0  0
 12 18  1  0  0  0  0
 12 19  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
5793

> <PUBCHEM_CONFORMER_RMSD>
0.6

> <PUBCHEM_CONFORMER_DIVERSEORDER>
1
8
4
3
6
2
7
5

> <PUBCHEM_MMFF94_PARTIAL_CHARGES>
17
1 -0.56
10 0.28
11 0.56
12 0.28
2 -0.68
20 0.4
21 0.4
22 0.4
23 0.4
24 0.4
3 -0.68
4 -0.68
5 -0.68
6 -0.68
7 0.28
8 0.28
9 0.28

> <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
2.2

> <PUBCHEM_PHARMACOPHORE_FEATURES>
12
1 1 acceptor
1 2 acceptor
1 2 donor
1 3 acceptor
1 3 donor
1 4 acceptor
1 4 donor
1 5 acceptor
1 5 donor
1 6 acceptor
1 6 donor
6 1 7 8 9 10 11 rings

> <PUBCHEM_HEAVY_ATOM_COUNT>
12

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
4

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
1

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
1

> <PUBCHEM_CONFORMER_ID>
000016A100000001

> <PUBCHEM_MMFF94_ENERGY>
27.5671

> <PUBCHEM_FEATURE_SELFOVERLAP>
60.963

> <PUBCHEM_SHAPE_FINGERPRINT>
12423570 1 11440868573817502106
16945 1 18338792415782022985
18185500 45 18410572903041433514
193761 8 18122905588311858698
21040471 1 18339079405607498497
23235685 24 18410853278337711372
2334 1 17834396004942031816
23402655 69 18267854013893658469
23552423 10 18120089751981225934
23559900 14 18270413797459679636
241688 4 17978794506319527824
2748010 2 18194117417474202524
5084963 1 17769667725561077874
5255222 1 18335976562245387725
528862 383 18261390018417527851
528886 8 18411416206831951673
53812653 166 18412823616180991792
63268167 104 18410859824078808769
66348 1 18411141294496697448

> <PUBCHEM_SHAPE_MULTIPOLES>
211.74
3.59
2.42
0.65
1.75
0.15
0
-0.75
-0.13
-0.7
0.1
-0.03
0.03
-0.1

> <PUBCHEM_SHAPE_SELFOVERLAP>
421.584

> <PUBCHEM_SHAPE_VOLUME>
123.3

> <PUBCHEM_COORDINATE_TYPE>
2
5
10

$$$$
"""
    
    MoleculeView(molFile: testData)
        .ignoresSafeArea(edges: .all)
    
}
