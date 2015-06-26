import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

public class PegJumping2 {

	private static final int MAX_TIME = 14500;
	private static final boolean DEBUG = false;
	private static final byte NONE = 0;

	private final long endTime = System.currentTimeMillis() + MAX_TIME;
	private int N, NN, N2;
	private State best = new State();
	private XorShift rnd = new XorShift();

	public String[] getMoves(int[] pegValue, String[] board) {
		{// input
			N = board.length;
			NN = N * N;
			N2 = N * 2;
		}
		boolean white[][] = new boolean[4][NN];
		boolean black[][] = new boolean[4][NN];
		for (int w = 0; w < 4; ++w) {
			for (int i = 0; i < N; ++i)
				white[w][i] = (i + w) % 2 == 0;
			for (int y = 1; y < N; ++y) {
				for (int x = 0; x < N; ++x) {
					white[w][getPos(y, x)] = !white[w][getPos(y - 1, x)];
				}
			}
			for (int y = (w <= 4 >> 1 ? 0 : 1); y < N; y += 2) {
				for (int x = 0; x < N; ++x) {
					black[w][getPos(y, x)] = !white[w][getPos(y, x)];
				}
			}
		}

		State init = new State(pegValue, board);
		TIME: for (int w = 0;; w = (w + 1) % 4) {
			List<State> states = new ArrayList<>();
			states.add(init);
			while (true) {
				List<State> next = new ArrayList<>();
				Collections.sort(states, (o1, o2) -> o2.score - o1.score);
				for (int i = 0, size = Math.min(states.size(), 10); i < size; ++i) {
					next.addAll(states.get(i).child(white[w], black[w]));
				}
				if (next.isEmpty())
					break;
				states = next;
			}
			for (State s : states) {
				s.getScore2(white[w], black[w]);
			}
			if (System.currentTimeMillis() >= endTime)
				break TIME;
		}

		{// output
			List<String> res = new ArrayList<>();
			while (true) {
				if (best.x == null)
					break;
				res.add(orderToString(best.n));
				best = best.x;
			}
			Collections.reverse(res);
			return res.toArray(new String[0]);
		}
	}

	private class State {
		State x;
		int[] n;
		byte[] s;
		int score;

		State() {
			score = 0;
		}

		State(int[] pegValue, String[] board) {
			x = null;
			n = null;
			s = new byte[NN];
			score = 0;

			for (int i = 0; i < NN; ++i) {
				char c = board[getY(i)].charAt(getX(i));
				if (c != '.')
					s[i] = (byte) pegValue[c - '0'];
			}
		}

		State(State x, int[] n, byte[] s, int score) {
			this.x = x;
			this.n = n;
			this.s = s;
			this.score = score;
		}

		private final boolean isStart(int p, byte[] s) {
			int x = getX(p);
			return (x + 2 < N && s[p + 1] != NONE && s[p + 2] == NONE)
					|| (x - 2 >= 0 && s[p - 1] != NONE && s[p - 2] == NONE)
					|| (p + N2 < NN && s[p + N] != NONE && s[p + N2] == NONE)
					|| (p - N2 >= 0 && s[p - N] != NONE && s[p - N2] == NONE);
		}

		private final int next(byte[] s, int p) {
			int res = 0, x = getX(p);
			if (x + 1 < N && s[p + 1] == NONE)
				res++;
			if (x - 1 >= 0 && s[p - 1] == NONE)
				res++;
			if (p + N < NN && s[p + N] == NONE)
				res++;
			if (p - N >= 0 && s[p - N] == NONE)
				res++;
			return res;
		}

		private final boolean isOK(byte[] s, int p, boolean ok[]) {
			int x = getX(p);
			return s[p] == NONE || (x + 2 < N && s[p + 1] != NONE && s[p + 2] == NONE && ok[p + 1] && ok[p + 2])
					|| (x - 2 >= 0 && s[p - 1] != NONE && s[p - 2] == NONE && ok[p - 1] && ok[p - 2])
					|| (p + N2 < NN && s[p + N] != NONE && s[p + N2] == NONE && ok[p + N] && ok[p + N2])
					|| (p - N2 >= 0 && s[p - N] != NONE && s[p - N2] == NONE && ok[p - N] && ok[p - N2]);
		}

		private final State move(State parent, byte[] s, int p, boolean ok[]) {
			int x = getX(p), np, dp;
			np = p + 2;
			dp = p + 1;
			if (x + 2 < N && s[dp] != NONE && s[np] == NONE && ok[dp] && ok[np]) {
				s[np] = s[p];
				s[p] = s[dp] = NONE;
				return new State(parent, new int[] { np, p }, s, score);
			}
			np = p - 2;
			dp = p - 1;
			if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE && ok[dp] && ok[np]) {
				s[np] = s[p];
				s[p] = s[dp] = NONE;
				return new State(parent, new int[] { np, p }, s, score);
			}
			np = p + N2;
			dp = p + N;
			if (np < NN && s[dp] != NONE && s[np] == NONE && ok[dp] && ok[np]) {
				s[np] = s[p];
				s[p] = s[dp] = NONE;
				return new State(parent, new int[] { np, p }, s, score);
			}
			np = p - N2;
			dp = p - N;
			if (np >= 0 && s[dp] != NONE && s[np] == NONE && ok[dp] && ok[np]) {
				s[np] = s[p];
				s[p] = s[dp] = NONE;
				return new State(parent, new int[] { np, p }, s, score);
			}
			throw new RuntimeException();
		}

		List<State> child(boolean[] white, boolean[] black) {
			List<State> res = new ArrayList<>();
			List<Integer> start = new ArrayList<>();
			for (int i = 0; i < NN; ++i) {
				if (white[i] && s[i] != NONE && isStart(i, s)) {
					start.add(i);
				}
			}
			if (start.isEmpty())
				return res;
			Set<Integer> used = new HashSet<>();
			for (int i = 0; i < 10; ++i) {
				int p = start.get(rnd.nextInt(start.size()));
				if (used.contains(p))
					continue;
				used.add(p);
				int x = getX(p), np, dp;
				np = p + 2;
				dp = p + 1;
				if (x + 2 < N && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p - 2;
				dp = p - 1;
				if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p + N2;
				dp = p + N;
				if (p + N2 < NN && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p - N2;
				dp = p - N;
				if (p - N2 >= 0 && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 16 : 0) + next(s, np)));
				}
			}
			return res;
		}

		void getScore2(boolean[] white, boolean[] black) {
			int next[] = getLongPath(s, white, black);
			byte[] s = Arrays.copyOf(this.s, this.s.length);
			boolean ok[] = new boolean[NN];
			Arrays.fill(ok, true);
			int score = 0;
			s[next[next.length - 1]] = s[next[0]];
			for (int k = 0; k < next.length - 1; ++k) {
				int a = next[k], b = next[k + 1];
				int deletePos = delPos(a, b);
				score += s[deletePos];
				s[deletePos] = NONE;
				ok[a] = ok[b] = ok[deletePos] = false;
			}
			boolean keep = true;
			State parent = this;
			while (keep) {
				keep = false;
				for (int k = 0; k < next.length; ++k) {
					{
						int p = next[k], px = getX(p), py = getY(p), tp1, tp2, tp3;
						for (int q = 0; q < 4; ++q) {
							if (q == 0) {
								if (px + 2 >= N || py + 2 >= N)
									continue;
								tp1 = p + 2;
								tp2 = p + 2 + N2;
								tp3 = p + N2;
							} else if (q == 1) {
								if (px - 2 < 0 || py + 2 >= N)
									continue;
								tp1 = p - 2;
								tp2 = p - 2 + N2;
								tp3 = p + N2;
							} else if (q == 2) {
								if (px + 2 >= N || py - 2 < 0)
									continue;
								tp1 = p + 2;
								tp2 = p + 2 - N2;
								tp3 = p - N2;
							} else {
								if (px - 2 < 0 || py - 2 < 0)
									continue;
								tp1 = p - 2;
								tp2 = p - 2 - N2;
								tp3 = p - N2;
							}
							int d1 = delPos(p, tp1), d2 = delPos(tp1, tp2), d3 = delPos(tp2, tp3), d4 = delPos(tp3, p);
							if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE && s[d4] != NONE) {
								ok[d1] = ok[d2] = ok[d3] = ok[d4] = false;
								if (!isOK(s, tp1, ok) || !isOK(s, tp2, ok) || !isOK(s, tp3, ok)) {
									ok[d1] = ok[d2] = ok[d3] = ok[d4] = true;
								} else {
									if (s[tp1] != NONE)
										parent = move(parent, s, tp1, ok);
									if (s[tp2] != NONE)
										parent = move(parent, s, tp2, ok);
									if (s[tp3] != NONE)
										parent = move(parent, s, tp3, ok);
									ok[tp1] = ok[tp2] = ok[tp3] = false;
									int tmp[] = Arrays.copyOf(next, next.length + 4);
									System.arraycopy(next, k, tmp, k + 4, next.length - k);
									tmp[k + 1] = tp1;
									tmp[k + 2] = tp2;
									tmp[k + 3] = tp3;
									tmp[k + 4] = p;
									next = tmp;
									score += s[d1] + s[d2] + s[d3] + s[d4];
									s[d1] = s[d2] = s[d3] = s[d4] = NONE;
									keep = true;
									break;
								}
							}
						}
					}
					if (k + 1 < next.length) {
						int p1 = next[k], p2 = next[k + 1], tp1, tp2;
						for (int q = 0; q < 2; ++q) {
							if (Math.abs(p1 - p2) == 2) {
								tp1 = p1 + (q == 0 ? N2 : -N2);
								tp2 = p2 + (q == 0 ? N2 : -N2);
								if ((q == 0 && tp1 >= NN) || (q != 0 && tp1 < 0))
									continue;
							} else {
								tp1 = p1 + (q == 0 ? 2 : -2);
								tp2 = p2 + (q == 0 ? 2 : -2);
								if ((q == 0 && getX(p1) + 2 >= N) || (q != 0 && getX(p1) - 2 < 0))
									continue;
							}
							int d1 = delPos(p1, tp1), d2 = delPos(tp1, tp2), d3 = delPos(tp2, p2);
							if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE) {
								ok[d1] = ok[d2] = ok[d3] = false;
								if (!isOK(s, tp1, ok) || !isOK(s, tp2, ok)) {
									ok[d1] = ok[d2] = ok[d3] = true;
								} else {
									if (s[tp1] != NONE)
										parent = move(parent, s, tp1, ok);
									if (s[tp2] != NONE)
										parent = move(parent, s, tp2, ok);
									ok[tp1] = ok[tp2] = false;
									int tmp[] = Arrays.copyOf(next, next.length + 2);
									System.arraycopy(next, k, tmp, k + 2, next.length - k);
									tmp[k + 1] = tp1;
									tmp[k + 2] = tp2;
									next = tmp;
									score += s[d1] + s[d2] + s[d3];
									s[d1] = s[d2] = s[d3] = NONE;
									s[(p1 + p2) >> 1] = 1;
									keep = true;
									break;
								}
							}
						}
					}
				}
			}
			score *= next.length - 1;
			if (best.score < score) {
				best = new State(parent, next, s, score);
			}
		}

		void getScore(boolean[] black) {
			final int max = NN >> 2, width = 40, mask = (1 << 6) - 1;
			int size[] = new int[max];
			int pos[][] = new int[max][width];
			int prev[][] = new int[max][width];
			long bit[] = new long[(NN >> 6) + 1];
			for (int i = 0; i < NN; ++i) {
				if (black[i] && s[i] != NONE && isStart(i, s)) {
					byte startCell = s[i];
					s[i] = NONE;
					Arrays.fill(size, 0);
					size[0] = 1;
					pos[0][0] = i;
					int maxj = 0;
					for (int j = 0, sj = 0; j < max - 1; ++j) {
						int wsize = Math.min(width, size[j]);
						boolean used[] = new boolean[wsize];
						for (int w = 0; sj + 4 <= width && w < wsize; ++w) {
							int k = rnd.nextInt(size[j]);
							while (used[k])
								k = (k + 1) % size[j];
							used[k] = true;
							if (!isStart(pos[j][k], s))
								continue;
							Arrays.fill(bit, 0);
							for (int a = j, pre = k, pre2 = prev[a][pre]; a > 0; pre = pre2, --a, pre2 = prev[a][pre]) {
								int x = delPos(pos[a][pre], pos[a - 1][pre2]);
								bit[x >> 6] |= 1L << (x & mask);
							}
							int np = pos[j][k], x = getX(np), delp, nextp;
							delp = np + 1;
							nextp = np + 2;
							if (x + 2 < N && s[delp] != NONE && s[nextp] == NONE
									&& (bit[delp >> 6] & (1L << (delp & mask))) == 0L) {
								pos[j + 1][sj] = nextp;
								prev[j + 1][sj] = k;
								++sj;
							}
							delp = np - 1;
							nextp = np - 2;
							if (x - 2 >= 0 && s[delp] != NONE && s[nextp] == NONE
									&& (bit[delp >> 6] & (1L << (delp & mask))) == 0L) {
								pos[j + 1][sj] = nextp;
								prev[j + 1][sj] = k;
								++sj;
							}
							delp = np + N;
							nextp = np + N2;
							if (nextp < NN && s[delp] != NONE && s[nextp] == NONE
									&& (bit[delp >> 6] & (1L << (delp & mask))) == 0L) {
								pos[j + 1][sj] = nextp;
								prev[j + 1][sj] = k;
								++sj;
							}
							delp = np - N;
							nextp = np - N2;
							if (nextp >= 0 && s[delp] != NONE && s[nextp] == NONE
									&& (bit[delp >> 6] & (1L << (delp & mask))) == 0L) {
								pos[j + 1][sj] = nextp;
								prev[j + 1][sj] = k;
								++sj;
							}
						}
						if (sj > 0) {
							size[j + 1] = sj;
							sj = 0;
							maxj = j + 1;
						} else
							break;
					}
					for (int j = 0; j < size[maxj]; ++j) {
						int next[] = new int[maxj + 1], ns = 0;
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						boolean ok[] = new boolean[NN];
						Arrays.fill(ok, true);
						int score = 0;
						for (int k = maxj, pre = j; k > 0; pre = prev[k][pre], --k) {
							int a = pos[k][pre], b = pos[k - 1][prev[k][pre]];
							next[ns++] = a;
							int deletePos = delPos(a, b);
							score += s[deletePos];
							s[deletePos] = NONE;
							ok[a] = ok[b] = ok[deletePos] = false;
						}
						next[ns++] = i;
						s[pos[maxj][j]] = startCell;
						boolean keep = true;
						State parent = this;
						while (keep) {
							keep = false;
							for (int k = 0; k < next.length; ++k) {
								{
									int p = next[k], px = getX(p), py = getY(p), tp1, tp2, tp3;
									for (int q = 0; q < 4; ++q) {
										if (q == 0) {
											if (px + 2 >= N || py + 2 >= N)
												continue;
											tp1 = p + 2;
											tp2 = p + 2 + N2;
											tp3 = p + N2;
										} else if (q == 1) {
											if (px - 2 < 0 || py + 2 >= N)
												continue;
											tp1 = p - 2;
											tp2 = p - 2 + N2;
											tp3 = p + N2;
										} else if (q == 2) {
											if (px + 2 >= N || py - 2 < 0)
												continue;
											tp1 = p + 2;
											tp2 = p + 2 - N2;
											tp3 = p - N2;
										} else {
											if (px - 2 < 0 || py - 2 < 0)
												continue;
											tp1 = p - 2;
											tp2 = p - 2 - N2;
											tp3 = p - N2;
										}
										int d1 = delPos(p, tp1), d2 = delPos(tp1, tp2), d3 = delPos(tp2, tp3), d4 = delPos(
												tp3, p);
										if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE && s[d4] != NONE) {
											ok[d1] = ok[d2] = ok[d3] = ok[d4] = false;
											if (!isOK(s, tp1, ok) || !isOK(s, tp2, ok) || !isOK(s, tp3, ok)) {
												ok[d1] = ok[d2] = ok[d3] = ok[d4] = true;
											} else {
												if (s[tp1] != NONE)
													parent = move(parent, s, tp1, ok);
												if (s[tp2] != NONE)
													parent = move(parent, s, tp2, ok);
												if (s[tp3] != NONE)
													parent = move(parent, s, tp3, ok);
												ok[tp1] = ok[tp2] = ok[tp3] = false;
												int tmp[] = Arrays.copyOf(next, next.length + 4);
												System.arraycopy(next, k, tmp, k + 4, next.length - k);
												tmp[k + 1] = tp1;
												tmp[k + 2] = tp2;
												tmp[k + 3] = tp3;
												tmp[k + 4] = p;
												next = tmp;
												score += s[d1] + s[d2] + s[d3] + s[d4];
												s[d1] = s[d2] = s[d3] = s[d4] = NONE;
												keep = true;
												break;
											}
										}
									}
								}
								if (k + 1 < next.length) {
									int p1 = next[k], p2 = next[k + 1], tp1, tp2;
									for (int q = 0; q < 2; ++q) {
										if (Math.abs(p1 - p2) == 2) {
											tp1 = p1 + (q == 0 ? N2 : -N2);
											tp2 = p2 + (q == 0 ? N2 : -N2);
											if ((q == 0 && tp1 >= NN) || (q != 0 && tp1 < 0))
												continue;
										} else {
											tp1 = p1 + (q == 0 ? 2 : -2);
											tp2 = p2 + (q == 0 ? 2 : -2);
											if ((q == 0 && getX(p1) + 2 >= N) || (q != 0 && getX(p1) - 2 < 0))
												continue;
										}
										int d1 = delPos(p1, tp1), d2 = delPos(tp1, tp2), d3 = delPos(tp2, p2);
										if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE) {
											ok[d1] = ok[d2] = ok[d3] = false;
											if (!isOK(s, tp1, ok) || !isOK(s, tp2, ok)) {
												ok[d1] = ok[d2] = ok[d3] = true;
											} else {
												if (s[tp1] != NONE)
													parent = move(parent, s, tp1, ok);
												if (s[tp2] != NONE)
													parent = move(parent, s, tp2, ok);
												ok[tp1] = ok[tp2] = false;
												int tmp[] = Arrays.copyOf(next, next.length + 2);
												System.arraycopy(next, k, tmp, k + 2, next.length - k);
												tmp[k + 1] = tp1;
												tmp[k + 2] = tp2;
												next = tmp;
												score += s[d1] + s[d2] + s[d3];
												s[d1] = s[d2] = s[d3] = NONE;
												s[(p1 + p2) >> 1] = 1;
												keep = true;
												break;
											}
										}
									}
								}
							}
						}
						score *= next.length - 1;
						if (best.score < score) {
							best = new State(parent, next, s, score);
						}
					}
					s[i] = startCell;
				}
			}
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

	private char getDir(int diff) {
		if (diff == 2)
			return 'R';
		else if (diff == N2)
			return 'D';
		else if (diff == -2)
			return 'L';
		// else if (diff == -N2)
		else
			return 'U';
	}

	private String orderToString(int[] order) {
		int pos = order[order.length - 1];
		StringBuilder s = new StringBuilder();
		s.append(getY(pos)).append(" ").append(getX(pos)).append(" ");
		for (int i = order.length - 1; i > 0; --i) {
			s.append(getDir(order[i - 1] - order[i]));
		}
		return s.toString();
	}

	private static final int delPos(int p1, int p2) {
		return (p1 + p2) >> 1;
	}

	private class MinCostFlow {
		class edge {
			int to, cap, cost, rev;

			edge(int to, int cap, int cost, int rev) {
				this.to = to;
				this.cap = cap;
				this.cost = cost;
				this.rev = rev;
			}
		}

		int V;
		ArrayList<edge> G[];
		int h[];
		int dist[];
		int prevv[], preve[];
		int go[];

		@SuppressWarnings("unchecked")
		public MinCostFlow(int V) {
			this.V = V;
			G = new ArrayList[V];
			for (int i = 0; i < V; i++)
				G[i] = new ArrayList<edge>();
			h = new int[V];
			dist = new int[V];
			prevv = new int[V];
			preve = new int[V];
			go = new int[V];
		}

		void add_edge(int from, int to, int cap, int cost) {
			G[from].add(new edge(to, cap, cost, G[to].size()));
			G[to].add(new edge(from, 0, -cost, G[from].size() - 1));
		}

		/*
		 * 負の辺に対応できない・・・
		 * よく考えよう
		 */
		int min_cost_flow(int s, int t, int f) {
			final int max = Integer.MAX_VALUE >> 1;
			class pair implements Comparable<pair> {
				int dist, v;

				pair(int dist, int v) {
					this.dist = dist;
					this.v = v;
				}

				@Override
				public int compareTo(pair paramT) {
					return dist - paramT.dist;
				}
			}
			int res = 0;
			Arrays.fill(h, 0);
			Arrays.fill(go, -1);
			while (f > 0) {
				Arrays.fill(dist, max);
				dist[s] = 0;
				PriorityQueue<pair> queue = new PriorityQueue<pair>();
				queue.add(new pair(0, s));
				while (!queue.isEmpty()) {
					pair p = queue.poll();
					int v = p.v;
					if (dist[v] < p.dist)
						continue;
					for (int i = 0; i < G[v].size(); i++) {
						edge e = G[v].get(i);
						if (e.cap == 0)
							continue;
						int d = dist[v] + e.cost + h[v] - h[e.to];
						if (dist[e.to] > d) {
							dist[e.to] = d;
							prevv[e.to] = v;
							preve[e.to] = i;
							queue.add(new pair(d, e.to));
						}
					}
				}
				if (dist[t] == max)
					return -1;
				for (int v = 0; v < V; v++) {
					h[v] += dist[v];
				}
				int d = f;
				for (int v = t; v != s; v = prevv[v]) {
					d = Math.min(d, G[prevv[v]].get(preve[v]).cap);
				}
				f -= d;
				res += d * h[t];
				for (int v = t; v != s; v = prevv[v]) {
					edge e = G[prevv[v]].get(preve[v]);
					go[prevv[v]] = v;
					e.cap -= d;
					G[v].get(e.rev).cap += d;
				}
			}
			return res;
		}
	}

	private final int[] getLongPath(final byte[] X, final boolean[] white, final boolean[] black) {
		List<Integer> wlist = new ArrayList<>();
		Map<Integer, Integer> to = new HashMap<>();
		for (int i = 0; i < NN; ++i) {
			if (X[i] != NONE && white[i]) {
				to.put(i, wlist.size());
				wlist.add(i);
			}
		}
		int V = 2 + wlist.size() * 2, s = V - 2, t = V - 1;
		MinCostFlow flow = new MinCostFlow(V);
		for (int i = 0; i < wlist.size(); ++i) {
			int p = wlist.get(i), x = getX(p), y = getY(p), outi = wlist.size() + 1;
			flow.add_edge(i, outi, 1, 0);
			flow.add_edge(outi, t, 1, 0);
			if (x - 1 >= 0 && x + 1 < N && black[p + 1]) {
				if ((X[p + 1] == NONE && X[p - 1] != NONE) || (X[p + 1] != NONE && X[p - 1] == NONE))
					flow.add_edge(s, i, 1, 0);
				if (X[p + 1] == NONE) {
					int p1, p2;
					p1 = p + 1 - N;
					p2 = p + 1 - N2;
					if (p2 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p + 1 + N;
					p2 = p + 1 + N2;
					if (p2 < NN && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p + 2;
					p2 = p + 3;
					if (x + 3 < N && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
				}
				if (X[p - 1] == NONE) {
					int p1, p2;
					p1 = p - 1 - N;
					p2 = p - 1 - N2;
					if (p2 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p - 1 + N;
					p2 = p - 1 + N2;
					if (p2 < NN && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p - 2;
					p2 = p - 3;
					if (x - 3 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
				}
			} else if (y - 1 >= 0 && y + 1 < N && black[p + N]) {
				if ((X[p + N] == NONE && X[p - N] != NONE) || (X[p + N] != NONE && X[p - N] == NONE))
					flow.add_edge(s, i, 1, 0);
				if (X[p + N] == NONE) {
					int p1, p2;
					p1 = p + N - 1;
					p2 = p + N - 2;
					if (x - 2 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p + N + 1;
					p2 = p + N + 2;
					if (x + 2 < N && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p + N + N;
					p2 = p + N + N2;
					if (p2 < NN && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
				}
				if (X[p - N] == NONE) {
					int p1, p2;
					p1 = p - N - 1;
					p2 = p - N - 2;
					if (x - 2 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p - N + 1;
					p2 = p - N + 2;
					if (x + 2 < N && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
					p1 = p - N - N;
					p2 = p - N - N2;
					if (p2 >= 0 && X[p1] != NONE && X[p2] == NONE) {
						flow.add_edge(outi, to.get(p1), 1, -1);
					}
				}
			}
		}
		System.out.println(1);
		flow.min_cost_flow(s, t, 1);
		List<Integer> tmp = new ArrayList<>();
		System.out.println(2);
		int now = flow.go[t];
		while (now != s) {
			if (now < wlist.size()) {
				tmp.add(wlist.get(now));
			}
			now = flow.go[now];
		}
		System.out.println(3);
		int res[] = new int[tmp.size() + 1];
		for (int i = 0; i < tmp.size() - 1; ++i) {
			int p1 = tmp.get(i), p2 = tmp.get(i + 1);
			if (p1 + N < NN && black[p1 + N]) {
				res[i] = p1 + (p1 < p2 ? N : -N);
			} else {
				res[i] = p1 + (p1 < p2 ? 1 : -1);
			}
		}
		res[0] = tmp.get(0) + tmp.get(0) - res[1];
		res[res.length - 1] = tmp.get(tmp.size() - 1) + tmp.get(tmp.size() - 1) + res[res.length - 2];
		return res;
	}
}
