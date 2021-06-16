package shoving;

import static java.lang.Math.sqrt;

public class Biomass {

	public Results pushing2D(int nBac, double[] bac_x, double[] bac_y, double[] bac_r, double tol) {

		int shov = 1; // boolean test for more shoving (1: shoving needed; 0: shoving ready)
		double[] shx = new double[nBac];
		double[] shy = new double[nBac];

		double dx, dy, d, overlap, dd, dxd, dyd;

		while (shov == 1) { // one shoving step
			shov = 0; // start assuming there is no need for pushing

			for (int i = 0; i < nBac; i++) { // loop over all n cells
				for (int j = i + 1; j < nBac; j++) { // search for overlap with any other cell j
					dx = bac_x[i] - bac_x[j];
					dy = bac_y[i] - bac_y[j];

					// calculate euclidean distance d between cell i (reference cell) and the new
					// cell j
					d = sqrt(dx * dx + dy * dy + 1e-20);

					// calculate overlap
					overlap = bac_r[i] + bac_r[j] - d;
					dd = overlap / d;
					dxd = dx * dd;
					dyd = dy * dd;
					if (dd > tol) {
						shx[i] += dxd;
						shy[i] += dyd;

						shx[j] -= dxd;
						shy[j] -= dyd;
						shov = 1;
					}
				}
			}

			for (int i = 0; i < nBac; i++) {
				// move cell i with the resultant components in x, y, z directions
				bac_x[i] = bac_x[i] + shx[i];
				bac_y[i] = bac_y[i] + shy[i];

				// reset the resultant vectors of movement to zero
				shx[i] = 0;
				shy[i] = 0;
			}
		}

		Results Results = new Results();
		Results.bac_x = bac_x;
		Results.bac_y = bac_y;
		return Results;
	}

	public Results pushing3D(int n, double[] bac_x, double[] bac_y, double[] bac_z, double[] bac_r) {

		int shov = 1; // boolean test for more shoving (1: shoving needed; 0: shoving ready)
		double[] shx = new double[n];
		double[] shy = new double[n];
		double[] shz = new double[n];

		double dx, dy, dz, d, overlap, dd, dxd, dyd, dzd;

		while (shov == 1) { // one shoving step
			shov = 0; // start assuming there is no need for pushing

			for (int i = 0; i < n; i++) { // loop over all n cells
				for (int j = i + 1; j < n; j++) { // search for overlap with any other cell j
					dx = bac_x[i] - bac_x[j];
					dy = bac_y[i] - bac_y[j];
					dz = bac_z[i] - bac_z[j];

					// calculate euclidian distance d between cell i (reference cell) and the new
					// cell j
					d = sqrt(dx * dx + dy * dy + dz * dz + 1e-20);

					// calculate overlap
					overlap = bac_r[i] + bac_r[j] - d;
					dd = overlap / d;
					dxd = dx * dd;
					dyd = dy * dd;
					dzd = dz * dd;
					if (dd > 0.1) {
						shx[i] = shx[i] + dxd;
						shy[i] = shy[i] + dyd;
						shz[i] = shz[i] + dzd;

						shx[j] = shx[j] - dxd;
						shy[j] = shy[j] - dyd;
						shz[j] = shz[j] - dzd;
						shov = 1;
					}
				}
			}

			for (int i = 0; i < n; i++) {
				// move cell i with the resultant components in x, y, z directions
				bac_x[i] = bac_x[i] + shx[i];
				bac_y[i] = bac_y[i] + shy[i];
				bac_z[i] = bac_z[i] + shz[i];

				// reset the resultant vectors of movement to zero
				shx[i] = 0;
				shy[i] = 0;
				shz[i] = 0;
			}
		}

		Results Results = new Results();
		Results.bac_x = bac_x;
		Results.bac_y = bac_y;
		Results.bac_z = bac_z;
		return Results;
	}

}
