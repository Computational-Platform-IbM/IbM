package shoving;

import java.util.List;

public class Bacterium {
	public double x;
	public double y;
	public double r;
	public double shx;
	public double shy;
	public int index;
	public List<Bacterium> neighbours;

	public Bacterium(int i, double x, double y, double r) {
		this.x = x;
		this.y = y;
		this.r = r;
		this.index = i;
		this.shx = 0;
		this.shy = 0;
	}
}
