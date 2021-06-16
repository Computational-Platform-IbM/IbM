package shoving;

public class Main {

	public static void main(String[] args) {

//		double w = 200;
//		double h = 200;
//		int nBac = 100;
//		double max_r = 5;
//
//		Rectangle canvas = new Rectangle(w / 2, h / 2, w / 2, h / 2);
//		Quadtree qt = new Quadtree(canvas, 4);
//
//		// generate bacteria
//		for (int i = 0; i < nBac; i++) {
//			double x = Math.random() * w;
//			double y = Math.random() * h;
//			double r = Math.random() * max_r;
//			Bacterium b = new Bacterium(i, x, y, r);
//			boolean success = qt.insert(b);
//			if (!success) {
//				System.out.printf("%d", i);
//			}
//		}
//
//		System.out.printf("%d bacteria in tree\n", qt.query(canvas, null).size());
//
//		Rectangle searchSpace = new Rectangle(20, 20, 20, 20);
//		List<Bacterium> results = qt.query(searchSpace, null);
//		System.out.printf("%d bacteria found in searchSpace\n", results.size());
//		for (Bacterium b : results) {
//			System.out.printf("%d: %.2f %.2f\n", b.index, b.x, b.y);
//
//		}

//		BiomassQuadtree bioqt = new BiomassQuadtree(0, 200, 0, 200);
//		Biomass bio = new Biomass();
//		double[] x_og = { 52.1196, 54.6674, 42.7789, 46.2535, 59.9449, 42.5632, 43.5799, 55.0585, 53.2432, 55.6862 };
//		double[] y_og = { 41.9379, 41.1714, 59.2479, 52.3311, 41.7326, 51.2254, 52.3305, 59.2769, 51.4861, 47.4232 };
//		double[] r_og = { 2.3564, 1.6056, 2.7079, 1.5853, 2.7511, 2.4289, 2.5534, 3.4693, 3.1967, 1.2072 };
//		double max_r = 5;
//		Results rQt = bioqt.pushing2D(10, x_og, y_og, r_og, 0.2, 4 * max_r);
//		System.out.printf("%.2f %.2f\n", rQt.bac_x[0], rQt.bac_y[0]);
//
//		double[] x_og2 = { 52.1196, 54.6674, 42.7789, 46.2535, 59.9449, 42.5632, 43.5799, 55.0585, 53.2432, 55.6862 };
//		double[] y_og2 = { 41.9379, 41.1714, 59.2479, 52.3311, 41.7326, 51.2254, 52.3305, 59.2769, 51.4861, 47.4232 };
//		double[] r_og2 = { 2.3564, 1.6056, 2.7079, 1.5853, 2.7511, 2.4289, 2.5534, 3.4693, 3.1967, 1.2072 };
//		Results rBio = bio.pushing2D(10, x_og2, y_og2, r_og2, 0.2);
//		System.out.printf("%.2f %.2f\n", rBio.bac_x[0], rBio.bac_y[0]);

//		System.out.printf("old position: %.2f %.2f\nnew position: %.2f %.2f\n", 52.1196, 41.9379, r.bac_x[0],
//				r.bac_y[0]);

	}
}
