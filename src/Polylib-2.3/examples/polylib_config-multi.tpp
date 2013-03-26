Polylib {
	sphere{
		filepath[@]="sphere.stl"
		
	}
	car {
		filepath = "./car.stl"
		movable = "true"
		velocity = 0.05
	}
	windmill {
		class_name  = "PolygonGroup"
		blades {
			movable = "true"
			center_x = 0.0
			center_y = 123.45
			center_z = 345.67
			filepath[@] = "./blade1.stl"
			filepath[@] = "./blade2.stl"
			filepath[@] = "./blade3.stl"
		}
		tower {
			class_name = "PolygonGroup"
			filepath = "./tower.stl"
		}
	}
} // end of Polylib