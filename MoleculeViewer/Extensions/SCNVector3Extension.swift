//
//  SCNVector3Extension.swift
//  MoleculeViewer
//
//  Created by Pierce Oxley on 23/1/25.
//

import Foundation
import SceneKit

extension SCNVector3 {
    func length() -> Float {
        return sqrt(x * x + y * y + z * z)
        
    }
    
    func normalized() -> SCNVector3 {
        let len = length()
        return SCNVector3(x / len, y / len, z / len)
    }
    
    func clamped(max: Float) -> SCNVector3 {
        let magnitude = self.length()
        return (magnitude > max) ? (self.normalized() * max) : self
    }
    
    func cross(_ v: SCNVector3) -> SCNVector3 {
        return SCNVector3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        )
    }
    
    func dot(_ v: SCNVector3) -> Float {
        return x * v.x + y * v.y + z * v.z
    }
    
    func perpendicular() -> SCNVector3 {
        // Return a vector perpendicular to this one (not unique)
        if x != 0 || y != 0 {
            return SCNVector3(-y, x, 0)
        } else {
            return SCNVector3(0, -z, y)
        }
    }
    
    static func + (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
        return SCNVector3(left.x + right.x, left.y + right.y, left.z + right.z)
    }
    
    static func - (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
        return SCNVector3(left.x - right.x, left.y - right.y, left.z - right.z)
    }
    
    static func * (vector: SCNVector3, scalar: Float) -> SCNVector3 {
        return SCNVector3(vector.x * scalar, vector.y * scalar, vector.z * scalar)
    }
    
    static func / (vector: SCNVector3, scalar: Float) -> SCNVector3 {
        return SCNVector3(vector.x / scalar, vector.y / scalar, vector.z / scalar)
    }
}
