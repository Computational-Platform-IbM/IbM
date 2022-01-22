package shoving;

import static java.lang.Math.sqrt;

import java.util.List;

public class BiomassQuadtree {

	private Rectangle boundary;

	public BiomassQuadtree(double minx, double maxx, double miny, double maxy) {
		double x = (maxx + minx) / 2;
		double y = (maxy + miny) / 2;
		double w = (maxx + minx) / 2;
		double h = (maxy + miny) / 2;
		this.boundary = new Rectangle(x, y, w, h);
	}

	public Results pushing2D(int nBac, double[] bac_x, double[] bac_y, double[] bac_r, double tol,
			double maxCheckingDistance, double kDist) {
		// create a quadtree
		Quadtree qt = new Quadtree(this.boundary, 10);

		Bacterium[] bacteria = new Bacterium[nBac];
		// first convert all x, y, r into bacterium instances and insert into quadtree
		for (int i = 0; i < nBac; i++) {
			bacteria[i] = new Bacterium(i, bac_x[i], bac_y[i], bac_r[i]);
			qt.insert(bacteria[i]);
		}

		// perform the shoving
		int shov = 1; // boolean test for more shoving (1: shoving needed; 0: shoving ready)
		double dx, dy, d, overlap, dd, dxd, dyd;

		// get neighbouring cells for each cell and store per cell (assumes that cells
		// will not move very far thus neighbouring cells stay constant).
		// To make sure that there are no edge cases, the search area is increased to
		// 4*max_radius as the quadtree is not updated after every shoving step (only
		// the bacteria)
		for (Bacterium b : bacteria) {
			Rectangle searchSpace = new Rectangle(b.x, b.y, maxCheckingDistance * 2, maxCheckingDistance * 2); // assume
			// maxCheckingDistance=2*max_radius
			List<Bacterium> neighbours = qt.query(searchSpace, null);
			b.neighbours = neighbours;
		}

		int cycle = 0;
		while (shov == 1) { // one shoving step
			cycle += 1;
			shov = 0; // start assuming there is no need for pushing

			for (Bacterium b1 : bacteria) { // loop over all n cells
				for (Bacterium b2 : b1.neighbours) { // only loop over cells in the neighbourhood of cell i
					if (b1 == b2 || b2.index < b1.index) { // make sure to not check cell with itself or to check
															// combinations twice
						continue;
					}
					dx = b1.x - b2.x;
					dy = b1.y - b2.y;

					// calculate euclidean distance d between cell i (reference cell) and the new
					// cell j
					d = sqrt(dx * dx + dy * dy + 1e-20);

					// calculate overlap
					overlap = kDist * (b1.r + b2.r) - d;
					dd = overlap / d;
					dxd = dx * dd;
					dyd = dy * dd;
					if (dd > tol) {
						b1.shx += dxd / 2;
						b1.shy += dyd / 2;

						b2.shx -= dxd / 2;
						b2.shy -= dyd / 2;
						shov = 1;
					}
				}
			}

			for (Bacterium b : bacteria) { // move bacteria and reset shx & shy
				// move cell i with the resultant components in x, y directions
				b.x += b.shx;
				b.y += b.shy;

				b.shx = 0;
				b.shy = 0;
			}
		}

		System.out.printf("\nshoving done in %d cycles...\n", cycle);

		// convert bacteria list into bac_x and bac_y again
		for (int i = 0; i < nBac; i++) {
			bac_x[i] = bacteria[i].x;
			bac_y[i] = bacteria[i].y;
		}

		Results Results = new Results();
		Results.bac_x = bac_x;
		Results.bac_y = bac_y;
		return Results;
	}

}
