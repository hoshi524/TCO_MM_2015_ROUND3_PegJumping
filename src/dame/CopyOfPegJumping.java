package dame;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CopyOfPegJumping {

	private static final int MAX_TIME = 14500;
	private final long endTime = System.currentTimeMillis() + MAX_TIME;
	private static final boolean DEBUG = true;
	private static final int NONE = -1;

	private int N, NN, N2;
	private State best;
	private char[] dirc = { 'R', 'D', 'L', 'U' };
	private XorShift rnd = new XorShift();

	private byte getDir(int diff) {
		if (diff == 2)
			return 0;
		if (diff == N2)
			return 1;
		if (diff == -2)
			return 2;
		if (diff == -N2)
			return 3;
		throw new RuntimeException();
	}

	public String[] getMoves(int[] pegValue, String[] board) {
		{// input
			N = board.length;
			NN = N * N;
			N2 = N * 2;
			best = new State();
		}
		boolean white[] = new boolean[NN];
		int w = 0;
		TIME: while (true) {
			List<State> states = new ArrayList<>();
			states.add(new State(pegValue, board));
			{// setWhite
				for (int i = 0; i < N; ++i)
					white[i] = (i + w) % 2 == 0;
				for (int y = 1; y < N; ++y) {
					for (int x = 0; x < N; ++x) {
						white[getPos(y, x)] = !white[getPos(y - 1, x)];
					}
				}
			}
			while (true) {
				List<State> next = new ArrayList<>();
				Collections.sort(states, (o1, o2) -> o2.score - o1.score);
				for (int i = 0, size = Math.min(states.size(), 15); i < size; ++i) {
					next.addAll(states.get(i).center(white));
				}
				if (next.isEmpty())
					break;
				states = next;
			}
			for (State s : states)
				s.getScore(white);
			if (System.currentTimeMillis() >= endTime)
				break TIME;
			++w;
			// break;
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

		State() {
			score = Integer.MIN_VALUE;
		}

		State(int[] pegValue, String[] board) {
			x = null;
			n = null;
			s = new byte[NN];
			score = 0;

			for (int i = 0; i < NN; ++i) {
				char c = board[getY(i)].charAt(getX(i));
				if (c == '.')
					s[i] = NONE;
				else
					s[i] = (byte) pegValue[c - '0'];
			}
		}

		State(State x, Next n, byte[] s, int score) {
			this.x = x;
			this.n = n;
			this.s = s;
			this.score = score;
		}

		List<State> center(boolean[] white) {
			List<State> res = new ArrayList<>();
			List<Integer> start = new ArrayList<>();
			for (int i = 0; i < NN; ++i) {
				if (white[i] && s[i] > NONE) {
					boolean ok = false;
					ok |= getX(i) + 2 < N && s[i + 1] != NONE && s[i + 2] == NONE;
					ok |= getX(i) - 2 >= 0 && s[i - 1] != NONE && s[i - 2] == NONE;
					ok |= i + N2 < NN && s[i + N] != NONE && s[i + N2] == NONE;
					ok |= i - N2 >= 0 && s[i - N] != NONE && s[i - N2] == NONE;
					if (ok)
						start.add(i);
				}
			}
			if (start.isEmpty())
				return res;

			Set<Integer> used = new HashSet<>();
			final int max = Math.min(1000, NN / 2), width = 30;
			int size[] = new int[max];
			int pos[][] = new int[max][width];
			byte dir[][] = new byte[max][width];
			byte prev[][] = new byte[max][width];
			for (int i = 0; i < 20; ++i) {
				{
					int p = start.get(rnd.nextInt(start.size()));
					if (used.contains(p))
						continue;
					used.add(p);
					Arrays.fill(size, 0);
					size[0] = 1;
					pos[0][0] = p;
				}
				int maxj = 0;
				byte[] s = Arrays.copyOf(this.s, this.s.length);
				for (int j = 0; j < max - 1; ++j) {
					for (byte k = 0; k < size[j]; ++k) {
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
					if (size[j + 1] > 0)
						maxj = j + 1;
				}
				if (maxj > 0) {
					Next next = new Next();
					next.pos = pos[0][0];
					s = Arrays.copyOf(this.s, this.s.length);
					for (int k = maxj, pre = 0; k > 0; pre = prev[k][pre], --k) {
						next.d.add(dir[k][pre]);
						int deletePos = (pos[k][pre] + pos[k - 1][prev[k][pre]]) / 2;
						s[deletePos] = NONE;
					}
					s[pos[maxj][0]] = s[pos[0][0]];
					s[pos[0][0]] = NONE;
					res.add(new State(this, next, s, this.score + maxj));
				}
			}
			return res;
		}

		void getScore(boolean[] white) {
			for (int i = 0; i < NN; ++i) {
				if (!white[i] && s[i] > NONE) {
					boolean ok = false;
					ok |= getX(i) + 2 < N && s[i + 1] != NONE && s[i + 2] == NONE;
					ok |= getX(i) - 2 >= 0 && s[i - 1] != NONE && s[i - 2] == NONE;
					ok |= i + N2 < NN && s[i + N] != NONE && s[i + N2] == NONE;
					ok |= i - N2 >= 0 && s[i - N] != NONE && s[i - N2] == NONE;
					if (ok) {
						List<Integer> path = new EulerPath().eulerPath(s, i);
						Next next = new Next();
						next.pos = i;
						int score = 0;
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						for (int k = 0, size = path.size() - 1; k < size - 1; ++k) {
							next.d.add(getDir(path.get(k + 1) - path.get(k)));
							int deletePos = (path.get(k) + path.get(k + 1)) / 2;
							score += s[deletePos];
							s[deletePos] = NONE;
						}
						s[path.get(0)] = s[i];
						s[i] = NONE;
						score *= path.size();
						if (best.score < score) {
							best = new State(this, next, s, this.score + score);
						}
					}
				}
			}
		}
	}

	private class EulerPath {
		void visit(byte[] s, int[][] adj, int p, List<Integer> path) {
			int x = getX(p);
			if (x + 2 < N && s[p + 1] != NONE && s[p + 2] == NONE && adj[p][p + 2] != 0) {
				--adj[p][p + 2];
				--adj[p + 2][p];
				visit(s, adj, p + 2, path);
			} else if (x - 2 >= 0 && s[p - 1] != NONE && s[p - 2] == NONE && adj[p][p - 2] != 0) {
				--adj[p][p - 2];
				--adj[p - 2][p];
				visit(s, adj, p - 2, path);

			} else if (p + N2 < NN && s[p + N] != NONE && s[p + N2] == NONE && adj[p][p + N2] != 0) {
				--adj[p][p + N2];
				--adj[p + N2][p];
				visit(s, adj, p + N2, path);
			} else if (p - N2 >= 0 && s[p - N] != NONE && s[p - N2] == NONE && adj[p][p - N2] != 0) {
				--adj[p][p - N2];
				--adj[p - N2][p];
				visit(s, adj, p - N2, path);
			}
			path.add(p);
		}

		List<Integer> eulerPath(byte[] s, int p) {
			byte tmp = s[p];
			s[p] = NONE;
			List<Integer> path = new ArrayList<>();
			int adj[][] = new int[NN][NN];
			for (int i = 0; i < s.length; ++i) {
				int x = getX(i);
				if (x + 2 < N && s[i + 1] != NONE && s[i + 2] == NONE)
					++adj[i][i + 2];
				if (x - 2 >= 0 && s[i - 1] != NONE && s[i - 2] == NONE)
					++adj[i][i - 2];
				if (i + N2 < NN && s[i + N] != NONE && s[i + N2] == NONE)
					++adj[i][i + N2];
				if (i - N2 >= 0 && s[i - N] != NONE && s[i - N2] == NONE)
					++adj[i][i - N2];
			}
			visit(s, adj, p, path);
			s[p] = tmp;
			return path;
		}
	}

	private class Next {
		int pos;
		List<Byte> d = new ArrayList<>();

		public String toString() {
			Collections.reverse(d);
			StringBuilder s = new StringBuilder();
			s.append(getY(pos)).append(" ").append(getX(pos)).append(" ");
			for (byte x : d)
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

	private final int getPos(int y, int x) {
		return y * N + x;
	}

	private static void debug(final Object... obj) {
		if (DEBUG)
			System.err.println(Arrays.deepToString(obj));
	}
}
