import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PegJumping {

	private static final int MAX_TIME = 14500;
	private final long endTime = System.currentTimeMillis() + MAX_TIME;
	private static final boolean DEBUG = true;
	private static final int NONE = -1;

	private int N, NN, N2;
	private State best = new State(new String[0]);
	private int[] dir, pegValue;
	private char[] dirc = { 'R', 'D', 'L', 'U' };
	private XorShift rnd = new XorShift();

	String[] getMoves(int[] pegValue, String[] board) {
		{// input
			N = board.length;
			NN = N * N;
			N2 = N * 2;
			dir = new int[] { 2, N2, -2, -N2 };
			this.pegValue = pegValue;
		}
		List<State> states = new ArrayList<>();
		states.add(new State(board));
		TIME: while (true) {
			while (true) {
				if (System.currentTimeMillis() >= endTime)
					break TIME;
				List<State> next = new ArrayList<>();
				for (int i = 0, size = Math.min(states.size(), 20); i < size; ++i) {
					State s = states.get(i);
					if (best.score < s.score)
						best = s;
					next.addAll(s.next());
				}
				if (next.isEmpty())
					break;
				states = next;
				Collections.sort(states, (o1, o2) -> o2.score - o1.score);
			}
			break;
		}

		{// output
			List<String> res = new ArrayList<>();
			while (true) {
				if (best.x == null)
					break;
				res.add(best.n.toString());
				best = best.x;
			}
			Collections.reverse(res);
			return res.toArray(new String[0]);
		}
	}

	private class State {
		State x;
		Next n;
		byte[] s;
		int score;

		State(String[] board) {
			x = null;
			n = null;
			s = new byte[NN];
			score = 0;

			for (int i = 0; i < NN; ++i) {
				char c = board[getY(i)].charAt(getX(i));
				if (c == '.')
					s[i] = NONE;
				else
					s[i] = (byte) (c - '0');
			}
		}

		State(State x, Next n, byte[] s, int score) {
			this.x = x;
			this.n = n;
			this.s = s;
			this.score = score;
		}

		List<State> next() {
			List<State> res = new ArrayList<>();
			List<Integer> start = new ArrayList<>();
			for (int i = 0; i < NN; ++i) {
				if (s[i] > NONE) {
					boolean ok = false;
					ok |= getX(i) + 2 < N && s[i + 1] != NONE && s[i + 2] == NONE;
					ok |= getX(i) - 2 >= 0 && s[i - 1] != NONE && s[i - 2] == NONE;
					ok |= i + N2 < NN && s[i + N] != NONE && s[i + N2] == NONE;
					ok |= i - N2 >= 0 && s[i - N] != NONE && s[i - N2] == NONE;
					if (ok) {
						start.add(i);
					}
				}
			}
			if (start.isEmpty())
				return res;
			Set<Integer> used = new HashSet<>();
			final int max = Math.min(1000, NN / 2), width = 40;
			int size[] = new int[max];
			int pos[][] = new int[max][width];
			int dir[][] = new int[max][width];
			int prev[][] = new int[max][width];
			for (int i = 0; i < 20; ++i) {
				{
					Arrays.fill(size, 0);
					int p = start.get(rnd.nextInt(start.size()));
					if (used.contains(p))
						continue;
					used.add(p);
					size[0] = 1;
					pos[0][0] = p;
				}
				byte[] s = Arrays.copyOf(this.s, this.s.length);
				for (int j = 0; j < max - 1; ++j) {
					for (int k = 0; k < size[j]; ++k) {
						int np = pos[j][k], x = getX(np);
						if (x + 2 < N && s[np + 1] != NONE && s[np + 2] == NONE) {
							pos[j + 1][size[j + 1]] = np + 2;
							dir[j + 1][size[j + 1]] = 0;
							prev[j + 1][size[j + 1]] = k;
							s[np + 1] = NONE;
							++size[j + 1];
						}
						if (x - 2 >= 0 && s[np - 1] != NONE && s[np - 2] == NONE) {
							pos[j + 1][size[j + 1]] = np - 2;
							dir[j + 1][size[j + 1]] = 2;
							prev[j + 1][size[j + 1]] = k;
							s[np - 1] = NONE;
							++size[j + 1];
						}
						if (np + N2 < NN && s[np + N] != NONE && s[np + N2] == NONE) {
							pos[j + 1][size[j + 1]] = np + N2;
							dir[j + 1][size[j + 1]] = 1;
							prev[j + 1][size[j + 1]] = k;
							s[np + N] = NONE;
							++size[j + 1];
						}
						if (np - N2 >= 0 && s[np - N] != NONE && s[np - N2] == NONE) {
							pos[j + 1][size[j + 1]] = np - N2;
							dir[j + 1][size[j + 1]] = 3;
							prev[j + 1][size[j + 1]] = k;
							s[np - N] = NONE;
							++size[j + 1];
						}
						if (size[j + 1] > width - 5)
							break;
					}
				}
				for (int j = max - 1; j > 0; --j) {
					if (size[j] > 0) {
						Next next = new Next();
						next.pos = pos[0][0];
						int score = 0;
						s = Arrays.copyOf(this.s, this.s.length);
						for (int k = j, pre = 0; k > 0; pre = prev[k][pre], --k) {
							next.d.add(dir[k][pre]);
							int deletePos = (pos[k][pre] + pos[k - 1][prev[k][pre]]) / 2;
							score += pegValue[s[deletePos]];
							s[deletePos] = NONE;
						}
						s[pos[j][0]] = s[pos[0][0]];
						s[pos[0][0]] = NONE;
						score *= j;
						Collections.reverse(next.d);
						res.add(new State(this, next, s, this.score + score));
						break;
					}
				}
			}
			return res;
		}
	}

	private class Next {
		int pos;
		List<Integer> d = new ArrayList<>();

		public String toString() {
			StringBuilder s = new StringBuilder();
			s.append(getY(pos)).append(" ").append(getX(pos)).append(" ");
			for (int x : d)
				s.append(dirc[x]);
			return s.toString();
		}
	}

	private final class XorShift {
		int x = 123456789;
		int y = 362436069;
		int z = 521288629;
		int w = 88675123;

		int nextInt(int n) {
			final int t = x ^ (x << 11);
			x = y;
			y = z;
			z = w;
			w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));
			final int r = w % n;
			return r < 0 ? r + n : r;
		}

		int nextInt() {
			final int t = x ^ (x << 11);
			x = y;
			y = z;
			z = w;
			w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));
			return w;
		}
	}

	private final int getY(int p) {
		return p / N;
	}

	private final int getX(int p) {
		return p % N;
	}

	private static void debug(final Object... obj) {
		if (DEBUG)
			System.err.println(Arrays.deepToString(obj));
	}
}
