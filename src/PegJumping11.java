import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PegJumping11 {

	private static final int MAX_TIME = 14500;
	private static final boolean DEBUG = false;
	private static final byte NONE = 0;

	private final long endTime = System.currentTimeMillis() + MAX_TIME;
	private int N, NN, N2;
	private State best = new State();
	private XorShift rnd = new XorShift();
	private long hash[][];

	public String[] getMoves(int[] pegValue, String[] board) {
		{// input
			N = board.length;
			NN = N * N;
			N2 = N * 2;
			hash = new long[NN][11];
			for (int i = 0; i < NN; ++i)
				for (int j = 0; j < 11; ++j)
					hash[i][j] = rnd.nextLong();
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
		int maxw = 0;
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
				State tmp = s.getScore(black[w]);
				if (best.score < tmp.score) {
					best = tmp;
					maxw = w;
				}
			}
			if (System.currentTimeMillis() >= endTime)
				break TIME;
		}
		//		while (true) {
		//			State tmp = best.getScore(black[maxw]);
		//			if (tmp.score == 0)
		//				break;
		//			best = tmp;
		//		}

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

		long hash() {
			long res = 0;
			for (int i = 0; i < NN; ++i)
				res ^= hash[i][s[i]];
			return res;
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

		State getScore(boolean[] black) {
			State res = new State();
			final int max = NN >> 2, width = 40, mask = (1 << 6) - 1;
			int size[] = new int[max], index[] = new int[NN], queue[] = new int[max], qs, qe;
			int dist[] = new int[NN], greedPrev[] = new int[NN];
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
						Arrays.fill(index, -1);
						for (int k = 0; k < next.length; ++k)
							index[next[k]] = k;
						while (keep) {
							keep = false;
							for (int k = 0; k < next.length; ++k) {
								int p = next[k];
								queue[0] = p;
								qs = 0;
								qe = 1;
								dist[p] = 0;
								Arrays.fill(greedPrev, -1);
								while (qs < qe) {
									int qp = queue[qs++];
									int x = getX(qp), np, dp;
									np = qp + 2;
									dp = qp + 1;
									if (x + 2 < N && s[dp] != NONE && s[np] == NONE && greedPrev[np] == -1) {
										greedPrev[np] = qp;
										dist[np] = dist[qp] + 1;
										if (index[np] == -1) {
											queue[qe++] = np;
										} else if (index[qp] == -1 && index[np] != -1
												&& Math.max(2, Math.abs(index[p] - index[np])) < dist[np]) {
											int add = dist[np] - Math.abs(index[p] - index[np]);
											int start = Math.min(index[p], index[np]);
											for (int w = start, end = Math.max(index[p], index[np]) - 1; w < end; ++w) {
												s[delPos(next[w], next[w + 1])] = 1;
												index[next[w + 1]] = -1;
											}
											int tmp[] = Arrays.copyOf(next, next.length + add);
											System.arraycopy(next, k, tmp, k + add, next.length - k);
											int tp = np, nindex = start + add + 1;
											while (tp != p) {
												int prevp = greedPrev[tp];
												s[delPos(tp, prevp)] = NONE;
												index[tp] = nindex;
												tmp[nindex] = tp;
												--nindex;
												tp = greedPrev[tp];
											}
											next = tmp;
											keep = true;
											break;
										}
									}
									np = qp - 2;
									dp = qp - 1;
									if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE && greedPrev[np] == -1) {
										greedPrev[np] = qp;
										dist[np] = dist[qp] + 1;
										if (index[np] == -1) {
											queue[qe++] = np;
										} else if (index[qp] == -1 && index[np] != -1
												&& Math.max(2, Math.abs(index[p] - index[np])) < dist[np]) {
											int add = dist[np] - Math.abs(index[p] - index[np]);
											int start = Math.min(index[p], index[np]);
											for (int w = start, end = Math.max(index[p], index[np]) - 1; w < end; ++w) {
												s[delPos(next[w], next[w + 1])] = 1;
												index[next[w + 1]] = -1;
											}
											int tmp[] = Arrays.copyOf(next, next.length + add);
											System.arraycopy(next, k, tmp, k + add, next.length - k);
											int tp = np, nindex = start + add + 1;
											while (tp != p) {
												int prevp = greedPrev[tp];
												s[delPos(tp, prevp)] = NONE;
												index[tp] = nindex;
												tmp[nindex] = tp;
												--nindex;
												tp = greedPrev[tp];
											}
											next = tmp;
											keep = true;
											break;
										}
									}
									np = qp + N2;
									dp = qp + N;
									if (np < NN && s[dp] != NONE && s[np] == NONE && greedPrev[np] == -1) {
										greedPrev[np] = qp;
										dist[np] = dist[qp] + 1;
										if (index[np] == -1) {
											queue[qe++] = np;
										} else if (index[qp] == -1 && index[np] != -1
												&& Math.max(2, Math.abs(index[p] - index[np])) < dist[np]) {
											int add = dist[np] - Math.abs(index[p] - index[np]);
											int start = Math.min(index[p], index[np]);
											for (int w = start, end = Math.max(index[p], index[np]) - 1; w < end; ++w) {
												s[delPos(next[w], next[w + 1])] = 1;
												index[next[w + 1]] = -1;
											}
											int tmp[] = Arrays.copyOf(next, next.length + add);
											System.arraycopy(next, k, tmp, k + add, next.length - k);
											int tp = np, nindex = start + add + 1;
											while (tp != p) {
												int prevp = greedPrev[tp];
												s[delPos(tp, prevp)] = NONE;
												index[tp] = nindex;
												tmp[nindex] = tp;
												--nindex;
												tp = greedPrev[tp];
											}
											next = tmp;
											keep = true;
											break;
										}
									}
									np = qp - N2;
									dp = qp - N;
									if (np >= 0 && s[dp] != NONE && s[np] == NONE && greedPrev[np] == -1) {
										greedPrev[np] = qp;
										dist[np] = dist[qp] + 1;
										if (index[np] == -1) {
											queue[qe++] = np;
										} else if (index[qp] == -1 && index[np] != -1
												&& Math.max(2, Math.abs(index[p] - index[np])) < dist[np]) {
											int add = dist[np] - Math.abs(index[p] - index[np]);
											int start = Math.min(index[p], index[np]);
											System.err.println(index[p] + " " + index[np]);
											for (int w = start, end = Math.max(index[p], index[np]) - 1; w < end; ++w) {
												s[delPos(next[w], next[w + 1])] = 1;
												index[next[w + 1]] = -1;
											}
											int tmp[] = Arrays.copyOf(next, next.length + add);
											System.arraycopy(next, k, tmp, k + add, next.length - k);
											int tp = np, nindex = start + add + 1;
											while (tp != p) {
												int prevp = greedPrev[tp];
												s[delPos(tp, prevp)] = NONE;
												index[tp] = nindex;
												tmp[nindex] = tp;
												--nindex;
												tp = greedPrev[tp];
											}
											next = tmp;
											keep = true;
											break;
										}
									}
								}
							}
						}
						score *= next.length - 1;
						if (res.score < score) {
							res = new State(parent, next, s, score);
						}
					}
					s[i] = startCell;
				}
			}
			return res;
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
			return w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));
		}

		long nextLong() {
			return ((long) nextInt() << 32) | (long) nextInt();
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
}
