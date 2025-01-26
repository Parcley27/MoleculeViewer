//
//  MoleculeView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SceneKit
import Foundation

struct Atom {
    var id: Int
    var position: SCNVector3
    var element: String
    var size: CGFloat
    var color: UIColor
    var bonds: [Bond] // Bonds directly associated with the atom
    
}

struct Bond {
    var atom1: Atom
    var atom2: Atom
    var type: Int // Bond type (e.g., single = 1, double = 2, etc.)
    var dipoleVector: SCNVector3? // Vector representing the dipole moment for this bond
    
}

struct MoleculeView: View {
    func loadAtomsAndBonds() -> [Atom]? {
        // Hardcoded .mol content including the header line
        let molContent = """
         24 25  0     0  0  0  0  0  0999 V2000
            0.4700    2.5688    0.0006 O   0  0  0  0  0  0  0  0  0  0  0  0
           -3.1271   -0.4436   -0.0003 O   0  0  0  0  0  0  0  0  0  0  0  0
           -0.9686   -1.3125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
            2.2182    0.1412   -0.0003 N   0  0  0  0  0  0  0  0  0  0  0  0
           -1.3477    1.0797   -0.0001 N   0  0  0  0  0  0  0  0  0  0  0  0
            1.4119   -1.9372    0.0002 N   0  0  0  0  0  0  0  0  0  0  0  0
            0.8579    0.2592   -0.0008 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.3897   -1.0264   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0307    1.4220   -0.0006 C   0  0  0  0  0  0  0  0  0  0  0  0
           -1.9061   -0.2495   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0
            2.5032   -1.1998    0.0003 C   0  0  0  0  0  0  0  0  0  0  0  0
           -1.4276   -2.6960    0.0008 C   0  0  0  0  0  0  0  0  0  0  0  0
            3.1926    1.2061    0.0003 C   0  0  0  0  0  0  0  0  0  0  0  0
           -2.2969    2.1881    0.0007 C   0  0  0  0  0  0  0  0  0  0  0  0
            3.5163   -1.5787    0.0008 H   0  0  0  0  0  0  0  0  0  0  0  0
           -1.0451   -3.1973   -0.8937 H   0  0  0  0  0  0  0  0  0  0  0  0
           -2.5186   -2.7596    0.0011 H   0  0  0  0  0  0  0  0  0  0  0  0
           -1.0447   -3.1963    0.8957 H   0  0  0  0  0  0  0  0  0  0  0  0
            4.1992    0.7801    0.0002 H   0  0  0  0  0  0  0  0  0  0  0  0
            3.0468    1.8092   -0.8992 H   0  0  0  0  0  0  0  0  0  0  0  0
            3.0466    1.8083    0.9004 H   0  0  0  0  0  0  0  0  0  0  0  0
           -1.8087    3.1651   -0.0003 H   0  0  0  0  0  0  0  0  0  0  0  0
           -2.9322    2.1027    0.8881 H   0  0  0  0  0  0  0  0  0  0  0  0
           -2.9346    2.1021   -0.8849 H   0  0  0  0  0  0  0  0  0  0  0  0
          1  9  2  0  0  0  0
          2 10  2  0  0  0  0
          3  8  1  0  0  0  0
          3 10  1  0  0  0  0
          3 12  1  0  0  0  0
          4  7  1  0  0  0  0
          4 11  1  0  0  0  0
          4 13  1  0  0  0  0
          5  9  1  0  0  0  0
          5 10  1  0  0  0  0
          5 14  1  0  0  0  0
          6  8  1  0  0  0  0
          6 11  2  0  0  0  0
          7  8  2  0  0  0  0
          7  9  1  0  0  0  0
         11 15  1  0  0  0  0
         12 16  1  0  0  0  0
         12 17  1  0  0  0  0
         12 18  1  0  0  0  0
         13 19  1  0  0  0  0
         13 20  1  0  0  0  0
         13 21  1  0  0  0  0
         14 22  1  0  0  0  0
         14 23  1  0  0  0  0
         14 24  1  0  0  0  0
        """
        
        let atomColors = ["H": UIColor(Color.white), "C": UIColor(Color.black), "N": UIColor(Color.blue), "O": UIColor(Color.red)]
        // From https://www.chem.ucla.edu/~harding/IGOC/A/atomic_radius.html
        let atomSizes = ["H": 0.25, "C": 0.70, "N": 0.65, "O": 0.60]
        
        let sizeScalar = 0.5
        
        
        var atoms: [Atom] = []
        var bonds: [Bond] = []
        
        let lines = molContent.split(separator: "\n")
        
        // Parse the header to determine the number of atoms and bonds
        guard let headerLine = lines.first else {
            print("Error: Missing .mol header.")
            
            return nil
            
        }
        
        let headerComponents = headerLine.split(separator: " ", omittingEmptySubsequences: true)
        
        guard headerComponents.count >= 2,
          let atomCount = Int(headerComponents[0]),
          let bondCount = Int(headerComponents[1]) else {
                print("Error: Invalid .mol header format.")
                
                return nil
            }
        
        // Parse atoms
        for i in 1...atomCount {
            let line = lines[i]
            let components = line.split(separator: " ", omittingEmptySubsequences: true)
            
            if components.count >= 4, let x = Float(components[0]), let y = Float(components[1]), let z = Float(components[2]) {
                let element = String(components[3])
                
                let atom = Atom(
                    id: i, // Atom ID corresponds to line index (1-based)
                    position: SCNVector3(x, y, z),
                    element: element,
                    size: (atomSizes[element] ?? 0.5) * sizeScalar,
                    color: atomColors[element] ?? UIColor(Color.gray),
                    bonds: []
                    
                )
                
                atoms.append(atom)
                
            }
        }
        
        // Parse bonds
        for i in (1 + atomCount)...(1 + atomCount + bondCount - 1) {
            let line = lines[i]
            let components = line.split(separator: " ", omittingEmptySubsequences: true)
            
            if components.count >= 3,
               let atom1Id = Int(components[0]),
               let atom2Id = Int(components[1]),
               let bondType = Int(components[2]) {
                
                let atom1 = atoms[atom1Id - 1] // 1-based index -> 0-based
                let atom2 = atoms[atom2Id - 1]
                
                // Create bond
                let bond = Bond(atom1: atom1, atom2: atom2, type: bondType)
                bonds.append(bond)
                
                // Update atom bonds
                atoms[atom1Id - 1].bonds.append(bond)
                atoms[atom2Id - 1].bonds.append(bond)
                
            }
        }
        
        return atoms
        
    }
    
    func createSphere(at position: SCNVector3, size: CGFloat = 0.1, color: UIColor = .gray) -> SCNNode {
        let sphere = SCNSphere(radius: size)
        
        let node = SCNNode(geometry: sphere)
        node.geometry?.firstMaterial?.diffuse.contents = color
        node.position = position
        
        return node
        
    }
    
    func createBonds(from atoms: [Atom]) -> [Bond] {
        var bonds: [Bond] = []
        
        // Set a threshold for bond distance (adjust based on the molecule)
        let bondThreshold: Float = 2.0
        
        for i in 0..<atoms.count {
            for j in i+1..<atoms.count {
                let atom1 = atoms[i]
                let atom2 = atoms[j]
                let distance = distanceBetween(atom1.position, atom2.position)
                
                if distance < bondThreshold {
                    // Create a bond with a default type of single bond
                    // Bond type could be updated later based on .mol data
                    let bond = Bond(atom1: atom1, atom2: atom2, type: 1)
                    bonds.append(bond)
                    
                }
            }
        }
        
        return bonds
    }
    
    func calculateDipoleVector(for bond: Bond) -> SCNVector3 {
        let electronegativities: [String: Float] = [
            "H": 2.20, "C": 2.55, "O": 3.44, "N": 3.04, "S": 2.58 // Example values
        ]
        
        // Get electronegativities of the two atoms
        let en1 = electronegativities[bond.atom1.element] ?? 0
        let en2 = electronegativities[bond.atom2.element] ?? 0
        
        // Calculate charge difference (q) and bond direction (d)
        let chargeDifference = abs(en1 - en2)
        let direction = SCNVector3(
            bond.atom2.position.x - bond.atom1.position.x,
            bond.atom2.position.y - bond.atom1.position.y,
            bond.atom2.position.z - bond.atom1.position.z
        )
        let distance = direction.length() // Bond length
        let normalizedDirection = direction.normalized() // Normalize the bond vector
        
        // Dipole vector = q * d
        let dipoleVector = normalizedDirection * chargeDifference * distance
        return dipoleVector
        
    }
    
    func calculateNetDipoleMoment(for bonds: [Bond]) -> SCNVector3 {
        var netDipoleMoment = SCNVector3(0, 0, 0)
        
        for bond in bonds {
            let dipoleVector = calculateDipoleVector(for: bond)
            netDipoleMoment = netDipoleMoment + dipoleVector
            
        }
        
        return netDipoleMoment
        
    }
    
    func createDipoleNode(from origin: SCNVector3 = SCNVector3(0, 0, 0), withDipole dipole: SCNVector3, color: UIColor = UIColor.magenta) -> SCNNode {
        let dipoleNode = SCNNode()
        
        // Calculate the length of the dipole vector
        let dipoleLength = dipole.length()
        
        // Create a cylinder to represent the dipole vector
        let cylinder = SCNCylinder(radius: 0.05, height: CGFloat(dipoleLength))
        let material = SCNMaterial()
        material.diffuse.contents = color // Set the arrow color
        cylinder.materials = [material]
        
        let cylinderNode = SCNNode(geometry: cylinder)
        
        // Position the cylinder at the midpoint of the dipole vector
        let midPoint = dipole * 0.5
        cylinderNode.position = midPoint
        
        // Align the cylinder with the dipole vector
        let up = SCNVector3(0, 1, 0) // Default up vector
        let axis = up.cross(dipole)  // Cross product for rotation axis
        let angle = acos(up.dot(dipole) / (up.length() * dipole.length())) // Angle between vectors
        
        if axis.length() > 0 {
            cylinderNode.rotation = SCNVector4(axis.x, axis.y, axis.z, angle)
        }
        
        // Add an arrowhead (optional)
        let cone = SCNCone(topRadius: 0, bottomRadius: 0.1, height: 0.2)
        let coneMaterial = SCNMaterial()
        coneMaterial.diffuse.contents = color
        cone.materials = [coneMaterial]
        
        let coneNode = SCNNode(geometry: cone)
        coneNode.position = dipole // Position the cone at the tip of the dipole
        coneNode.rotation = cylinderNode.rotation // Match the cylinder's rotation
        
        // Add the cylinder and cone to the dipole node
        dipoleNode.addChildNode(cylinderNode)
        dipoleNode.addChildNode(coneNode)
        
        return dipoleNode
    }
    
    func calculateDipoleDipoleInteraction(dipole1: SCNVector3, dipole2: SCNVector3, position1: SCNVector3, position2: SCNVector3) -> Float {
        // Distance vector between the two dipoles
        let rVector = SCNVector3(
            position2.x - position1.x,
            position2.y - position1.y,
            position2.z - position1.z
        )
        
        let r = rVector.length() // Magnitude of the distance vector
        let rUnit = rVector.normalized() // Unit vector in the direction of r
        
        // Dot products
        let mu1DotMu2 = dipole1.dot(dipole2)
        let mu1DotR = dipole1.dot(rUnit)
        let mu2DotR = dipole2.dot(rUnit)
        
        // Calculate the interaction energy
        let energy = -mu1DotMu2 / pow(r, 3) + 3 * mu1DotR * mu2DotR / pow(r, 3)
        
        return energy
    }
    
    func distanceBetween(_ point1: SCNVector3, _ point2: SCNVector3) -> Float {
        let dx = point2.x - point1.x
        let dy = point2.y - point1.y
        let dz = point2.z - point1.z
        
        return sqrt(dx * dx + dy * dy + dz * dz)
        
    }
    
    func createBondNode(from position1: SCNVector3, to position2: SCNVector3, bondType: Int = 1, color: UIColor = UIColor.gray, radius: CGFloat = 0.04) -> SCNNode {
        let bondNode = SCNNode()
        
        // Calculate the direction vector (from position1 to position2)
        let direction = SCNVector3(position2.x - position1.x, position2.y - position1.y, position2.z - position1.z)
        let midPoint = SCNVector3(
            (position1.x + position2.x) / 2,
            (position1.y + position2.y) / 2,
            (position1.z + position2.z) / 2
        )
        
        // Normalize the direction vector
        let normalizedDirection = direction.normalized()
                
        for bondNumber in 0..<bondType {
            // Calculate offset position
            //let bondOffset = offset * Float((bondNumber + bondType - 1)) // Spread bonds symmetrically
            //let bondOffset = offset * (Float((bondNumber + bondType)) - 2.5)
            
            // Perpendicular offset vector for multiple bonds
            var offset = normalizedDirection.perpendicular().normalized() * (0.2)
            
            if bondType == 3 {
                offset = normalizedDirection.perpendicular().normalized() * (0.15)
                
            }
            
            var bondOffset: SCNVector3 = .init()
            
            if bondType == 1 {
                bondOffset = offset * Float((bondNumber + bondType - 1))
                
            } else if bondType == 2 {
                bondOffset = offset * (Float((bondNumber + bondType)) - 2.5)
                
            } else {
                bondOffset = offset * Float((bondNumber + bondType - 4))
                
            }
            
            
            // Create a cylinder for the bond
            let length = distanceBetween(position1, position2)
            let cylinder = SCNCylinder(radius: radius, height: CGFloat(length))
            
            // Apply color to the cylinder
            let material = SCNMaterial()
            material.diffuse.contents = color
            cylinder.materials = [material]
            
            // Create a node for the cylinder
            let cylinderNode = SCNNode(geometry: cylinder)
            
            // Set position of the bond at the midpoint plus the offset
            cylinderNode.position = midPoint + bondOffset
            
            // Compute the rotation axis (cross product of the cylinder's local up vector and the direction vector)
            let up = SCNVector3(0, 1, 0)
            let axis = up.cross(direction)
            let angle = acos(up.dot(direction) / (up.length() * direction.length()))
            
            // Apply the rotation to align the cylinder
            cylinderNode.rotation = SCNVector4(axis.x, axis.y, axis.z, angle)
            
            bondNode.addChildNode(cylinderNode)
        }
        
        return bondNode
    }
    
    func createCylinder(from pointA: SCNVector3, to pointB: SCNVector3, color: UIColor = UIColor.gray, radius: CGFloat = 0.025) -> SCNNode {
        // Calculate the vector between pointA and pointB
        let direction = SCNVector3(pointB.x - pointA.x, pointB.y - pointA.y, pointB.z - pointA.z)
        
        // Calculate the distance between the two points (cylinder height)
        let height = distanceBetween(pointA, pointB)
        
        // Create the cylinder geometry
        let cylinder = SCNCylinder(radius: radius, height: CGFloat(height))
        
        // Apply color to the cylinder's material
        let material = SCNMaterial()
        material.diffuse.contents = color
        cylinder.materials = [material]
        
        // Create a node with the cylinder geometry
        let cylinderNode = SCNNode(geometry: cylinder)
        
        // Position the cylinder at the midpoint of pointA and pointB
        let midPoint = SCNVector3(
            (pointA.x + pointB.x) / 2,
            (pointA.y + pointB.y) / 2,
            (pointA.z + pointB.z) / 2
        )
        cylinderNode.position = midPoint
        
        // Align the cylinder along the vector from pointA to pointB
        let up = SCNVector3(0, 1, 0) // Default cylinder up direction
        let axis = up.cross(direction) // Rotation axis (cross product)
        let angle = acos(up.dot(direction) / (up.length() * direction.length())) // Rotation angle
        
        if axis.length() > 0 { // Avoid invalid rotations
            cylinderNode.rotation = SCNVector4(axis.x, axis.y, axis.z, angle)
        }
        
        return cylinderNode
        
    }
    
    func createMolecularStructure(at origin: SCNVector3 = SCNVector3(0, 0, 0), rotation: SCNVector3 = SCNVector3(0, 0, 0)) -> (node: SCNNode, dipole: SCNVector3, position: SCNVector3) {
        let structureNode = SCNNode() // Parent node for the molecular structure
        
        // Load atoms and bonds
        guard let atoms = loadAtomsAndBonds() else {
            print("Failed to load atoms and bonds")
            return (structureNode, SCNVector3(0, 0, 0), origin)
        }
        
        var allBonds: [Bond] = [] // Collect all bonds for dipole moment calculation
        
        // Create and add atom nodes
        for atom in atoms {
            let atomNode = createSphere(at: atom.position, size: atom.size, color: atom.color)
            structureNode.addChildNode(atomNode) // Add atom as child of structureNode
        }
        
        // Create and add bond nodes
        for atom in atoms {
            for bond in atom.bonds {
                // Only create bonds where this atom is the "first" to avoid duplication
                if bond.atom1.id < bond.atom2.id {
                    let bondNode = createBondNode(from: bond.atom1.position, to: bond.atom2.position, bondType: bond.type)
                    structureNode.addChildNode(bondNode) // Add bond as child of structureNode
                    
                    // Store the bond for dipole moment calculation
                    allBonds.append(bond)
                }
            }
        }
        
        // Calculate the net dipole moment
        let netDipoleMoment = calculateNetDipoleMoment(for: allBonds)
        print("Net Dipole Moment: \(netDipoleMoment)")
        print("Dipome Magnitude: \(sqrt(pow(netDipoleMoment.x, 2) + pow(netDipoleMoment.y, 2) + pow(netDipoleMoment.z, 2)))")

        
        
        // Add the dipole visualization
        let dipoleNode = createDipoleNode(from: origin, withDipole: netDipoleMoment, color: .blue)
        structureNode.addChildNode(dipoleNode)
        
        // Position and rotate the structure
        structureNode.position = origin
        structureNode.eulerAngles = SCNVector3(rotation.x, rotation.y, rotation.z)
        
        return (structureNode, netDipoleMoment, origin)
        
    }
    
    func createScene() -> SCNScene {
        let scene = SCNScene()
        
        /*
        let molecularStructureNode = createMolecularStructure(
            at: SCNVector3(0, 0, 0),
            //rotation: SCNVector3(0, Double.pi / 4, 0) // Rotate 45 degrees around the Y-axis
            rotation: SCNVector3(0, 0, 0)
            
        )
         */
        
        //scene.rootNode.addChildNode(molecularStructureNode)
        //scene.rootNode.addChildNode(createMolecularStructure())
        //scene.rootNode.addChildNode(createMolecularStructure(at: SCNVector3(0, 0, 10)))
        
        /*
        SCNTransaction.begin()
        SCNTransaction.animationDuration = 2.0 // Rotate over 2 seconds
         
        // Apply rotation
        molecularStructureNode.eulerAngles = SCNVector3(0, Double.pi, 0) // Rotate 180° around Y-axis
         
        SCNTransaction.commit()
        */
        
        let molecule1 = createMolecularStructure(at: SCNVector3(0, 5, 0), rotation: SCNVector3(0, 0, 0))
        let molecule2 = createMolecularStructure(at: SCNVector3(0, -5, 0), rotation: SCNVector3(0, 0, 0))
        
        scene.rootNode.addChildNode(molecule1.node)
        scene.rootNode.addChildNode(molecule2.node)
        
        // Calculate dipole-dipole interaction energy
        let interactionEnergy = calculateDipoleDipoleInteraction(
            dipole1: molecule1.dipole,
            dipole2: molecule2.dipole,
            position1: molecule1.position,
            position2: molecule2.position
        )
        
        print("Dipole-Dipole Interaction Energy: \(interactionEnergy)")
        
        // XYZ Coordniate reference
        let cylinderNodeX = createCylinder(from: SCNVector3(0, 0, 0), to: SCNVector3(1, 0, 0), color: .red)
        let cylinderNodeY = createCylinder(from: SCNVector3(0, 0, 0), to: SCNVector3(0, 1, 0), color: .green)
        let cylinderNodeZ = createCylinder(from: SCNVector3(0, 0, 0), to: SCNVector3(0, 0, 1), color: .blue)
        
        scene.rootNode.addChildNode(cylinderNodeX)
        scene.rootNode.addChildNode(cylinderNodeY)
        scene.rootNode.addChildNode(cylinderNodeZ)
        
        let lightNode = SCNNode()
        let light = SCNLight()
        light.type = .directional
        light.intensity = 1000
        light.color = UIColor.white
        light.castsShadow = true
        
        lightNode.light = light
        lightNode.eulerAngles = SCNVector3(-1, 0, Double.pi/4)
        lightNode.position = SCNVector3(0, 10, 0)
        scene.rootNode.addChildNode(lightNode)
        
        let ambientLightNode = SCNNode()
        let ambientLight = SCNLight()
        ambientLight.type = .ambient
        ambientLight.intensity = 200
        ambientLight.color = UIColor.white
        
        ambientLightNode.light = ambientLight
        scene.rootNode.addChildNode(ambientLightNode)
        
        let cameraNode = SCNNode()
        let camera = SCNCamera()
        camera.fieldOfView = 60
        camera.zNear = 1
        camera.zFar = 100
        camera.wantsDepthOfField = false
        camera.focalLength = 50
        //camera.aperture = 0.5

        cameraNode.camera = camera
        cameraNode.position = SCNVector3(0, 10, 50)
        cameraNode.look(at: SCNVector3(0, 0, 0))
        scene.rootNode.addChildNode(cameraNode)
        
        return scene
        
    }
    
    var body: some View {
        SceneView(
            scene: createScene(),
            options: [.allowsCameraControl],
            preferredFramesPerSecond: 120
        )
        .ignoresSafeArea()
        
    }
}

#Preview {
    MoleculeView()
    
}

