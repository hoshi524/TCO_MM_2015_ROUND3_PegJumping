import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PegJumping4 {

	private static final int MAX_TIME = 14500;
	private final long endTime = System.currentTimeMillis() + MAX_TIME;
	private static final boolean DEBUG = false;
	private static final byte NONE = 0;

	private int N, NN, N2;
	private State best;
	boolean white[][];
	boolean black[][];
	State[] type = new State[4];

	public String[] getMoves(int[] pegValue, String[] board) {
		{// input
			N = board.length;
			NN = N * N;
			N2 = N * 2;
			best = new State();
		}
		white = new boolean[4][NN];
		black = new boolean[4][NN];
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
		for (int i = 0; i < 4; ++i)
			type[i] = new State(pegValue, board);
		while (true) {
			for (int i = 0; i < 4; ++i) {
				get(i).getScore(black[i]);
			}
			if (System.currentTimeMillis() >= endTime)
				break;
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

	State get(int i) {
		if (type[i] == null)
			return null;
		State s = type[i];
		while (true) {
			List<State> child = s.child(white[i], black[i]);
			if (child.isEmpty()) {
				State res = s;
				type[i] = s.x;
				return res;
			} else if (s.ci >= child.size()) {
				s = s.x;
				if (s == null) {
					type[i] = null;
					return null;
				}
			} else {
				s = child.get(s.ci++);
			}
		}
	}

	private class State {
		State x;
		int[] n;
		byte[] s;
		int score, ci = 0;
		List<State> child = null;

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

		List<State> child(boolean[] white, boolean[] black) {
			if (child != null)
				return child;
			List<State> res = new ArrayList<>();
			for (int p = 0; p < NN; ++p) {
				if (white[p] && s[p] != NONE) {
					int x = getX(p), np, dp;
					np = p + 2;
					dp = p + 1;
					if (x + 2 < N && s[dp] != NONE && s[np] == NONE) {
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						s[np] = s[p];
						s[p] = s[dp] = NONE;
						res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 1 : 0)));
					}
					np = p - 2;
					dp = p - 1;
					if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE) {
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						s[np] = s[p];
						s[p] = s[dp] = NONE;
						res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 1 : 0)));
					}
					np = p + N2;
					dp = p + N;
					if (p + N2 < NN && s[dp] != NONE && s[np] == NONE) {
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						s[np] = s[p];
						s[p] = s[dp] = NONE;
						res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 1 : 0)));
					}
					np = p - N2;
					dp = p - N;
					if (p - N2 >= 0 && s[dp] != NONE && s[np] == NONE) {
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						s[np] = s[p];
						s[p] = s[dp] = NONE;
						res.add(new State(this, new int[] { np, p }, s, this.score + (black[dp] ? 1 : 0)));
					}
				}
			}
			Collections.sort(res, (o1, o2) -> o2.score - o1.score);
			return child = res;
		}

		void getScore(boolean[] black) {
			final int max = NN >> 2, width = 20, mask = (1 << 6) - 1;
			int size[] = new int[max];
			int pos[][] = new int[max][width];
			byte prev[][] = new byte[max][width];
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
						for (byte k = 0; sj + 4 <= width && k < size[j]; ++k) {
							if (!isStart(pos[j][k], s))
								continue;
							Arrays.fill(bit, 0);
							for (int a = j, pre = k, pre2 = prev[a][pre]; a > 0; pre = pre2, --a, pre2 = prev[a][pre]) {
								int x = (pos[a][pre] + pos[a - 1][pre2]) >> 1;
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
						int score = 0;
						for (int k = maxj, pre = j; k > 0; pre = prev[k][pre], --k) {
							int a = pos[k][pre], b = pos[k - 1][prev[k][pre]];
							next[ns++] = a;
							int deletePos = (a + b) >> 1;
							score += s[deletePos];
							s[deletePos] = NONE;
						}
						next[ns++] = i;
						s[pos[maxj][j]] = startCell;
						boolean keep = true;
						while (keep) {
							keep = false;
							// 貪欲に延ばす
							// まだまだパターンありそう
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
										if (s[tp1] != NONE || s[tp2] != NONE || s[tp3] != NONE)
											continue;
										int d1 = (p + tp1) >> 1;
										int d2 = (tp1 + tp2) >> 1;
										int d3 = (tp2 + tp3) >> 1;
										int d4 = (tp3 + p) >> 1;
										if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE && s[d4] != NONE) {
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
										if (s[tp1] != NONE || s[tp2] != NONE)
											continue;
										int d1 = (p1 + tp1) >> 1;
										int d2 = (tp1 + tp2) >> 1;
										int d3 = (tp2 + p2) >> 1;
										if (s[d1] != NONE && s[d2] != NONE && s[d3] != NONE) {
											int tmp[] = Arrays.copyOf(next, next.length + 2);
											System.arraycopy(next, k, tmp, k + 2, next.length - k);
											tmp[k + 1] = tp1;
											tmp[k + 2] = tp2;
											next = tmp;
											s[d1] = s[d2] = s[d3] = NONE;
											s[(p1 + p2) >> 1] = 1;
											keep = true;
											break;
										}
									}
								}
							}
						}
						score *= next.length - 1;
						if (best.score < this.score + score)
							best = new State(this, next, s, this.score + score);
					}
					s[i] = startCell;
				}
			}
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
}
