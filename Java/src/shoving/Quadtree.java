package shoving;

import java.util.ArrayList;
import java.util.List;

public class Quadtree {
	private Rectangle boundary;
	private int capacity;
	private Bacterium[] points;
	private int inserted;
	private Quadtree topleft;
	private Quadtree topright;
	private Quadtree bottomleft;
	private Quadtree bottomright;
	private boolean divided;

	public Quadtree(Rectangle boundary, int capacity) {
		this.boundary = boundary;
		this.capacity = capacity;
		this.points = new Bacterium[capacity];
		this.inserted = 0;
		this.divided = false;
	}

	public boolean insert(Bacterium b) {
		// Function to insert a bacterium into the quadtree. Returns true is successful.
		if (!this.boundary.contains(b)) {
			return false;
		}

		if (this.inserted < this.capacity) {
			this.points[this.inserted] = b;
			this.inserted += 1;
			return true;
		} else {
			if (!(this.divided)) {
				this.subdivide();
				this.divided = true;
			}
			// try to insert into every quadrant, if not possible, will continue to the next
			// (conditional method chaining)
			return (this.topleft.insert(b) || this.topright.insert(b) || this.bottomleft.insert(b)
					|| this.bottomright.insert(b));

		}
	}

	public List<Bacterium> query(Rectangle r, List<Bacterium> found) {
		if (found == null) {
			found = new ArrayList<Bacterium>();
		}

		// Function to return each of the bacteria in a certain region.
		if (!this.boundary.intersects(r)) {
			return found;
		} else {
			for (int j = 0; j < this.inserted; j++) {
				Bacterium p = points[j];
				if (r.contains(p)) {
					found.add(p);
				}
			}
			if (this.divided) {
				this.topleft.query(r, found);
				this.topright.query(r, found);
				this.bottomleft.query(r, found);
				this.bottomright.query(r, found);
			}
		}
		return found;
	}

	private void subdivide() {
		// Function to create new subtrees if the 'mother-tree' is full.
		Rectangle topleft = new Rectangle(this.boundary.x - this.boundary.w / 2, this.boundary.y + this.boundary.h / 2,
				this.boundary.w / 2, this.boundary.h / 2);
		Rectangle topright = new Rectangle(this.boundary.x + this.boundary.w / 2, this.boundary.y + this.boundary.h / 2,
				this.boundary.w / 2, this.boundary.h / 2);
		Rectangle bottomleft = new Rectangle(this.boundary.x - this.boundary.w / 2,
				this.boundary.y - this.boundary.h / 2, this.boundary.w / 2, this.boundary.h / 2);
		Rectangle bottomright = new Rectangle(this.boundary.x + this.boundary.w / 2,
				this.boundary.y - this.boundary.h / 2, this.boundary.w / 2, this.boundary.h / 2);
		this.topleft = new Quadtree(topleft, this.capacity);
		this.topright = new Quadtree(topright, this.capacity);
		this.bottomleft = new Quadtree(bottomleft, this.capacity);
		this.bottomright = new Quadtree(bottomright, this.capacity);
	}

	public Bacterium[] viewPoints() {
		return points;
	}
}
