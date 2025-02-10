//
//  MoleculeSimulationView.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 21/1/25.
//

import SwiftUI
import SceneKit
import Foundation

struct MoleculeSimulationView: View {
    let spawnArea: Int = 10
    let spawnSpeed: Double = 1.0
    let moleculeCount: Int = 30
    
    let sizeScalar: Double = 0.5
    
    let temperature: Float = 0.1
    let rotationRandomness: Float = Float.pi / 200
    let movementLoss: Float = 0.0001
    
    let originPull: Float = 10000
    let imfScalar: Float = 2
    
    func loadAtomsAndBonds() -> [Atom]? {
        // Hardcoded .mol content including the header line
        let molContent = """
        2  1  0     0  0  0  0  0  0999 V2000
         -0.5560    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
          0.5560    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        1  2  3  0  0  0  0
        """
        
        // From https://en.wikipedia.org/wiki/CPK_coloring
        let atomColors = ["H": UIColor(Color.white), "C": UIColor(Color.black), "N": UIColor(Color.blue), "O": UIColor(Color.red), "Na": #colorLiteral(red: 0.2745098174, green: 0.4862745106, blue: 0.1411764771, alpha: 1) , "Cl": UIColor(Color.green)]
        
        // From https://www.chem.ucla.edu/~harding/IGOC/A/atomic_radius.html
        let atomSizes = ["H": 0.25, "C": 0.70, "N": 0.65, "O": 0.60, "Na": 1.80, "Cl": 1.00]
        
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
                    size: (atomSizes[element] ?? 50) * sizeScalar,
                    color: atomColors[element] ?? UIColor(Color.pink),
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
    
    func randomJitter(scale: Float) -> SCNVector3 {
        let x = (Float.random(in: -1...1) * scale)
        let y = (Float.random(in: -1...1) * scale)
        let z = (Float.random(in: -1...1) * scale)
        return SCNVector3(x, y, z)
        
    }
    
    // MARK: Dipoles
    func calculateDipoleVector(for bond: Bond) -> SCNVector3 {
            
        // From https://www.chem.ucla.edu/~harding/IGOC/E/electronegativity.html
        let electronegativities: [String: Float] = [
            "H": 2.1, "C": 2.5, "N": 3.0, "O": 3.5, "Na": 0.9, "S": 2.5, "Cl": 3.0
            
        ]
        
        guard let en1 = electronegativities[bond.atom1.element],
              let en2 = electronegativities[bond.atom2.element] else {
            print("Warning: Electronegativity value missing for one or both atoms.")
            return SCNVector3(0, 0, 0)
            
        }
        
        let deltaEN = abs(en1 - en2)
        
        // Approximate partial charge (q) in elementary charge (e)
        let elementaryCharge: Float = 1.6e-19 // Coulombs
        let partialCharge = (deltaEN / 4.40) * elementaryCharge
        
        // Compute bond direction vector
        let direction = SCNVector3(
            bond.atom2.position.x - bond.atom1.position.x,
            bond.atom2.position.y - bond.atom1.position.y,
            bond.atom2.position.z - bond.atom1.position.z
            
        )
        
        let bondLength = direction.length() // Bond length in angstroms (assumption)
        let normalizedDirection = direction.normalized()
        
        // Convert bond length to meters (1 angstrom = 1e-10 m)
        let bondLengthMeters = bondLength * 1e-10
        
        // Compute dipole moment (C·m)
        let dipoleMoment = partialCharge * bondLengthMeters
        
        // Convert to Debye (1 D = 3.336e-30 C·m)
        let dipoleMomentDebye = dipoleMoment / 3.336e-30 // Debye
        print("Dipole moment: \(dipoleMomentDebye)")
        
        // Scale the vector based on dipole moment in Debye
        let dipoleVector = normalizedDirection * dipoleMomentDebye * -1
        
        return dipoleVector
        
    }
    
    func calculateNetDipoleMoment(for bonds: [Bond]) -> SCNVector3 {
        var netDipoleMoment = SCNVector3(0, 0, 0)
        
        let bondSlice = bonds.count / 2
        
        let bondSet = Array(bonds.prefix(bondSlice))
        
        for bond in bondSet {
            let dipoleVector = calculateDipoleVector(for: bond)
            print("Vector: \(dipoleVector)")
            netDipoleMoment = netDipoleMoment + dipoleVector
            
        }
        
        print("Net vector: \(netDipoleMoment) of \(netDipoleMoment.length()) D")
        
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
        
        // Physical constant (vacuum permittivity)
        let epsilon0: Float = 8.854e-12
        
        // Calculate the interaction energy with correct r^5 term
        let energy = (1 / (4 * Float.pi * epsilon0)) * ((-mu1DotMu2 / pow(r, 3)) + (3 * mu1DotR * mu2DotR / pow(r, 5)))
        
        return energy
        
    }
    
    func calculateAllDipoleInteractions(molecules: [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)]) -> [[Float]] {
        var interactionMatrix: [[Float]] = Array(repeating: Array(repeating: 0.0, count: molecules.count), count: molecules.count)
        
        for i in 0..<molecules.count {
            for j in i+1..<molecules.count {
                let interactionEnergy = calculateDipoleDipoleInteraction(
                    dipole1: molecules[i].dipole,
                    dipole2: molecules[j].dipole,
                    position1: molecules[i].position,
                    position2: molecules[j].position
                )
                interactionMatrix[i][j] = interactionEnergy
                interactionMatrix[j][i] = interactionEnergy // Symmetric matrix
                
            }
        }
        
        return interactionMatrix
        
    }
    
    // MARK: Forces
    func calculateNetForces(molecules: [(dipole: SCNVector3, position: SCNVector3)], repulsionConstant: Float = 12) -> [SCNVector3] {
        var forces = Array(repeating: SCNVector3(0, 0, 0), count: molecules.count)
        
        for i in 0..<molecules.count {
            let position1 = molecules[i].position
            let dipole1 = molecules[i].dipole
            
            // Attractive force towards origin
            let distanceFromOrigin = position1.length()
            if distanceFromOrigin > 0 { // Avoid division by zero
                let attractionMagnitude = pow(distanceFromOrigin, 2) / originPull
                let attractionForce = position1.normalized() * -attractionMagnitude
                forces[i] = forces[i] + attractionForce
                
            }
            
            for j in 0..<molecules.count where i != j {
                let dipole2 = molecules[j].dipole
                let position2 = molecules[j].position
                
                // Distance vector and magnitude
                let rVector = SCNVector3(
                    position2.x - position1.x,
                    position2.y - position1.y,
                    position2.z - position1.z
                )
                
                let r = rVector.length()
                let rUnit = rVector.normalized()
                
                // Dipole-dipole interaction force components
                let mu1DotMu2 = dipole1.dot(dipole2)
                let mu1DotR = dipole1.dot(rUnit)
                let mu2DotR = dipole2.dot(rUnit)
                
                let commonTerm = 3 * mu1DotR * mu2DotR - mu1DotMu2
                let dipoleForceMagnitude = commonTerm / pow(r, 4) // Gradient of dipole interaction potential
                
                // Repulsion force magnitude
                let repulsionForceMagnitude = repulsionConstant / pow(r, 6)
                
                // Total force magnitude (dipole + repulsion)
                let totalForceMagnitude = (dipoleForceMagnitude - repulsionForceMagnitude) * imfScalar
                
                // Force vector
                let force = rUnit * totalForceMagnitude
                forces[i] = forces[i] + force
                
            }
        }
        
        return forces
        
    }
    
    //MARK: Vector Display
    func visualizeNetForceVectors(
        scene: SCNScene,
        molecules: [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)]
    ) {
        // Remove old force vectors if re-visualizing
        let oldForceNodes = scene.rootNode.childNodes.filter { $0.name == "forceVector" }
        oldForceNodes.forEach { $0.removeFromParentNode() }
        
        // Calculate forces and visualize
        let forces = calculateNetForces(molecules: molecules.map { ($0.dipole, $0.position) })
        
        for i in 0..<molecules.count {
            let force = forces[i]
            let origin = molecules[i].position
            
            //print(molecules[i].position)
            let forceNode = createCylinder(from: origin, to: origin + force, color: .yellow)
            
            //let forceNode = createDipoleNode(withDipole: force, color: .yellow)
            forceNode.name = "forceVector" // Tag for easy removal later
            scene.rootNode.addChildNode(forceNode)
            
            //print("Molecule \(i + 1) Force: \(force)")
            
        }
    }
    
    func visualizeSingleForceVector(
        scene: SCNScene,
        bond: Bond
    ) {
        
        let forceVector = calculateDipoleVector(for: bond)
        
        let position1 = bond.atom1.position
        let position2 = bond.atom2.position
            
        let bondOrigin = SCNVector3(
                (position1.x + position2.x) / 2,
                (position1.y + position2.y) / 2,
                (position1.z + position2.z) / 2
            )
        
        let dipoleNode = createDipoleNode(from: bondOrigin, withDipole: forceVector, color: .red)
        scene.rootNode.addChildNode(dipoleNode)
        
    }
    
    func calculateCenterOfMass(for atoms: [Atom]) -> SCNVector3 {
        var totalMass: Float = 0
        var weightedSum = SCNVector3(0, 0, 0)
        
        let atomicMasses = ["H": 1, "C": 12, "N": 14, "O": 16, "Na": 23 , "Cl": 35]
        
        for atom in atoms {
            let mass = Float(atomicMasses[atom.element] ?? 1)
            totalMass += mass
            weightedSum = weightedSum + (atom.position * mass)
            
        }
        
        return totalMass > 0 ? weightedSum / totalMass : SCNVector3(0, 0, 0)
        
    }
    
    // MARK: Rotation
    func pointMolecules(
        scene: SCNScene,
        molecules: [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)]
    ) {
        for i in 0..<molecules.count {
            let molecule = molecules[i]
            let dipoleField = calculateLocalDipoleField(for: i, molecules: molecules)
            orientMolecule(
                node: molecule.node,
                localDipole: molecule.dipole,
                direction: dipoleField,
                rotationRandomness: rotationRandomness
                
            )
        }
    }
    
    func calculateLocalDipoleField(
        for moleculeIndex: Int,
        molecules: [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)]
    ) -> SCNVector3 {
        var field = SCNVector3Zero
        let referenceMolecule = molecules[moleculeIndex]
        
        for (index, otherMolecule) in molecules.enumerated() where index != moleculeIndex {
            // Compute world dipole of other molecule
            let worldDipole = rotateVector(otherMolecule.dipole, by: otherMolecule.node.rotation)
            let rVector = otherMolecule.position - referenceMolecule.position
            let distance = rVector.length()
            
            if distance > 0 {
                field = field + (worldDipole / (distance * distance)) // Inverse square scaling
                
            }
        }
        
        return field.normalized()
        
    }
    
    func orientMolecule(
        node: SCNNode,
        localDipole: SCNVector3,
        direction: SCNVector3,
        maxRotationAngle: CGFloat = .pi,
        rotationRandomness: Float = .pi / 180
    ) {
        guard localDipole.length() > 0 else { return }
        
        if localDipole.length() < 2 { return }
        
        // Get current world dipole direction
        let currentWorldDipole = rotateVector(localDipole, by: node.rotation).normalized()
        let desiredDirection = direction.normalized()
        
        guard currentWorldDipole.length() > 0 && desiredDirection.length() > 0 else { return }
        
        print(currentWorldDipole)
        
        // Calculate rotation axis and angle
        let axis = currentWorldDipole.cross(desiredDirection)
        guard axis.length() > 0 else { return }
        
        var angle = acos(currentWorldDipole.dot(desiredDirection))
        angle = angle.isNaN ? 0 : angle
        
        print(angle)
        
        // Apply jitter and clamp rotation
        //angle += Float.random(in: -rotationRandomness...rotationRandomness)
        let clampedAngle = min(angle, Float(maxRotationAngle))
        
        print(axis)
        print()
        
        // Apply incremental rotation
        node.rotation = SCNVector4(axis.x, axis.y, axis.z, clampedAngle)
        
    }
    
    func rotateVector(_ vector: SCNVector3, by rotation: SCNVector4) -> SCNVector3 {
        let quaternion = simd_quatf(angle: rotation.w, axis: simd_float3(rotation.x, rotation.y, rotation.z))
        let simdVector = simd_float3(vector.x, vector.y, vector.z)
        let rotated = quaternion.act(simdVector)
        return SCNVector3(rotated.x, rotated.y, rotated.z)
    }
    
    func updateMoleculePositions(
        molecules: inout [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)],
        deltaTime: Float
    ) {
        let forces = calculateNetForces(molecules: molecules.map { ($0.dipole, $0.position) })
        
        for i in 0..<molecules.count {
            let mass = molecules[i].mass
            let acceleration = forces[i] * (1.0 / mass) // F = ma => a = F/m
            
            molecules[i].velocity = (molecules[i].velocity * (1 - movementLoss)) + acceleration * deltaTime
            molecules[i].position = molecules[i].position + molecules[i].velocity * deltaTime
            
        }
    }
    
    func calculateTorque(forces: [SCNVector3], atoms: [Atom], centerOfMass: SCNVector3) -> SCNVector3 {
        atoms.enumerated().reduce(SCNVector3Zero) { torque, indexAtom in
            let (index, atom) = indexAtom
            let r = atom.position - centerOfMass // Relative to CoM
            return torque + r.cross(forces[index])
        }
    }
    
    // MARK: simulateMolecularMotion
   func simulateMolecularMotion(scene: SCNScene, molecules: inout [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)]) {
       let deltaTime: Float = 1.0 / 60.0 // 60 FPS
       
       // Create a local copy of the molecules to work with
       var moleculesCopy = molecules
        
       visualizeNetForceVectors(scene: scene, molecules: molecules)
       
       pointMolecules(scene: scene, molecules: molecules)
       
       //updateMoleculeRotationWithTorque(molecules: &molecules, deltaTime: deltaTime)
       
       // Update molecule positions and velocities
       updateMoleculePositions(molecules: &moleculesCopy, deltaTime: deltaTime)
       
       for i in 0..<molecules.count {
           let molecule = molecules[i]
           
           //print("Molecule \(i) position: \(molecule.node.position) rotation: \(molecule.node.rotation)")
           
           // Check if the position is too far or invalid
           if molecule.node.position.length().isNaN || molecule.node.position.length() > 1000 {
               print("Molecule \(i) in invalid position")
               exit(1)
               
           }
           
           // Check if the rotation contains NaN
           if molecule.node.rotation.x.isNaN || molecule.node.rotation.y.isNaN || molecule.node.rotation.z.isNaN {
               print("Molecule \(i) has invalid rotation!")
               
           }
       }
       
       // Animate molecule movement
       for i in 0..<moleculesCopy.count {
           let jitter = randomJitter(scale: temperature * 0.05) // Scale jitter by temperature
           molecules[i].position = molecules[i].position + jitter
           molecules[i].node.position = molecules[i].position
           
           let molecule = moleculesCopy[i]
           
           // Update node's position based on the molecule's new position
           let displacement = SCNVector3(
               x: molecule.position.x - molecule.node.position.x,
               y: molecule.position.y - molecule.node.position.y,
               z: molecule.position.z - molecule.node.position.z
           )
           
           // Create an SCNAction to move the molecule over time
           let moveAction = SCNAction.move(by: displacement, duration: TimeInterval(deltaTime))
           
           // Apply the move action to the molecule's node
           molecule.node.runAction(moveAction)
           
           // Update the position after applying the action
           moleculesCopy[i].node.position = molecule.position
       }
       
       // Add termination condition (e.g., stop when molecules move out of bounds)
       let maxDistance: Float = 1000.0
       if moleculesCopy.contains(where: { $0.position.length() > maxDistance }) {
           print("Simulation stopped: molecules moved out of bounds.")
           stopAnimationLoop()
           
       }
       
       // Reassign the modified copy back to the original inout parameter when needed
       molecules = moleculesCopy
       
   }
    
    @State private var molecules: [(node: SCNNode, dipole: SCNVector3, position: SCNVector3, velocity: SCNVector3, mass: Float)] = []
    @State private var timer: Timer?
    
    // Function to start the animation loop
    func startAnimationLoop(scene: SCNScene) {
        // Set up the timer for animation loop (60 FPS)
        timer = Timer.scheduledTimer(withTimeInterval: 1.0 / 10.0, repeats: true) { _ in
            // Update the simulation
            simulateMolecularMotion(scene: scene, molecules: &molecules)
            
        }
    }
    
    // Function to stop the animation loop
    func stopAnimationLoop() {
        timer?.invalidate()
        timer = nil
        
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
   
   // MARK: createMolecularStructure
   func createMolecularStructure(at origin: SCNVector3 = SCNVector3(0, 0, 0), rotation: SCNVector3 = SCNVector3(0, 0, 0)) -> (node: SCNNode, dipole: SCNVector3, position: SCNVector3) {
       let structureNode = SCNNode() // Parent node for the molecular structure
       
       // Load atoms and bonds
       guard var atoms = loadAtomsAndBonds() else {
           print("Failed to load atoms and bonds")
           return (structureNode, SCNVector3(0, 0, 0), origin)
       }
       
       let centerOfMass = calculateCenterOfMass(for: atoms)
       for i in 0..<atoms.count {
           atoms[i].position = atoms[i].position - centerOfMass
           
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
                   let bondNode = createBondNode(from: bond.atom1.position - centerOfMass, to: bond.atom2.position - centerOfMass, bondType: bond.type)
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
   
   // MARK: Scene
   func createScene() -> SCNScene {
       let scene = SCNScene()
       
       let molecule1Position = SCNVector3(-1.5, 0, 0)
       let molecule2Position = SCNVector3(1.5, 0, 0)
       
       var molecules = [
           (node: createMolecularStructure(at: molecule1Position, rotation: SCNVector3(0, 0, 0)).node,
            dipole: SCNVector3(x: -1.3493353, y: -0.9989537, z: 0.705799),
            position: molecule1Position,
            velocity: SCNVector3(0.1, 0, 0),
            mass: Float(0.1)),
            
           (node: createMolecularStructure(at: molecule2Position, rotation: SCNVector3(0, 0, 0)).node,
            dipole: SCNVector3(x: -1.3493353, y: -0.9989537, z: 0.705799),
            position: molecule2Position,
            velocity: SCNVector3(-0.1, 0, 0),
            mass: Float(0.1)),
           
        ]
       
        molecules.removeAll()
       
        for _ in 1 ... moleculeCount {
            let moleculeX = Int.random(in: -spawnArea ... spawnArea)
            let moleculeY = Int.random(in: -spawnArea ... spawnArea)
            let moleculeZ = Int.random(in: -spawnArea ... spawnArea)
           
            let moleculePosition = SCNVector3(moleculeX, moleculeY, moleculeZ)
           
            let speedX = Double.random(in: -spawnSpeed ... spawnSpeed)
            let speedY = Double.random(in: -spawnSpeed ... spawnSpeed)
            let speexZ = Double.random(in: -spawnSpeed ... spawnSpeed)
           
            let moleculeSpeed = SCNVector3(speedX, speedY, speexZ)
            
            let moleculeXRotation = Double.random(in: -.pi ... .pi)
            let moleculeYRotation = Double.random(in: -.pi ... .pi)
            let moleculeZRotation = Double.random(in: -.pi ... .pi)
            
            let moleculeOrientation = SCNVector3(moleculeXRotation, moleculeYRotation, moleculeZRotation)
           
            molecules.append(
            (node: createMolecularStructure(at: moleculePosition, rotation: moleculeOrientation).node,
             dipole: SCNVector3(x: -1.3493353, y: -0.9989537, z: 0.705799),
             position: moleculePosition,
             velocity: moleculeSpeed,
             mass: Float(0.1)))
           
        }
       
        // Add molecules to the scene
        for molecule in molecules {
            scene.rootNode.addChildNode(molecule.node)
           
        }
       
       // Start simulation
       simulateMolecularMotion(scene: scene, molecules: &molecules)
       
       timer = Timer.scheduledTimer(withTimeInterval: 1.0 / 60.0, repeats: true) { _ in
//           //Update the simulation
           simulateMolecularMotion(scene: scene, molecules: &molecules)
//
       }
       
        let _ = calculateAllDipoleInteractions(molecules: molecules)
        //print("Dipole-Dipole Interaction Matrix: \(interactionMatrix)")
       
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
        cameraNode.position = SCNVector3(0, 0, 30)
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
        .onAppear {
            // Start the animation loop when the view appears
            startAnimationLoop(scene: createScene())
           
        }
        .onDisappear {
            // Stop the animation loop when the view disappears
            stopAnimationLoop()
           
        }
    }
}

// MARK: Preview
#Preview {
    MoleculeSimulationView()
   
}
