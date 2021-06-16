package shoving;

public class Rectangle {
	public double x;
	public double y;
	public double w; // width defined as distance between left boundary and x-coordinate at center
						// (i.e. radius of the square)
	public double h; // height similarly defined

	public Rectangle(double x, double y, double w, double h) {
		this.x = x;
		this.y = y;
		this.w = w;
		this.h = h;
	}

	public boolean contains(Bacterium b) {
		// Function that returns true if the bacterium is present in this rectangle.
		return (b.x > this.x - this.w && b.x <= this.x + this.w && b.y > this.y - this.h && b.y <= this.y + this.h);
	}

	public boolean intersects(Rectangle r) {
		// Function that returns true if the given rectangle intersects this rectangle.
		return !(r.x - r.w > this.x + this.w || r.x + r.w < this.x - this.w || r.y - r.h > this.y + this.h
				|| r.y + r.h < this.y - this.h);
	}
}
