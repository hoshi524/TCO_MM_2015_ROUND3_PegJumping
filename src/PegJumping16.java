import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PegJumping16 {

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
				for (int i = 0, size = Math.min(states.size(), 2); i < size; ++i) {
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

		// best.debug();
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

		void debug() {
			for (int y = 0; y < N; ++y) {
				for (int x = 0; x < N; ++x) {
					System.out.print(String.format("%2d", s[getPos(y, x)]));
				}
				System.out.println();
			}
			System.out.println();
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
			int min = Math.min(10, start.size());
			boolean used[] = new boolean[start.size()];
			for (int i = 0; i < min; ++i) {
				int ri = rnd.nextInt(start.size());
				if (used[ri])
					ri = (ri + 1) % start.size();
				used[ri] = true;
				int p = start.get(ri), x = getX(p), np, dp;
				np = p + 2;
				dp = p + 1;
				if (x + 2 < N && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p - 2;
				dp = p - 1;
				if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p + N2;
				dp = p + N;
				if (p + N2 < NN && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, score + (black[dp] ? 16 : 0) + next(s, np)));
				}
				np = p - N2;
				dp = p - N;
				if (p - N2 >= 0 && s[dp] != NONE && s[np] == NONE) {
					byte[] s = Arrays.copyOf(this.s, this.s.length);
					s[np] = s[p];
					s[p] = s[dp] = NONE;
					res.add(new State(this, new int[] { np, p }, s, score + (black[dp] ? 16 : 0) + next(s, np)));
				}
			}
			return res;
		}

		State getScore(boolean[] black) {
			State res = new State();
			final int max = NN >> 1;
			int dist[] = new int[NN], index[] = new int[NN], queue[] = new int[max], prev[] = new int[NN];
			for (int i = 0; i < NN; ++i) {
				if (black[i] && s[i] != NONE && isStart(i, s)) {
					for (int q = 0; q < 20; ++q) {
						byte[] s = Arrays.copyOf(this.s, this.s.length);
						int next[];
						byte startCell = s[i];
						{
							List<Integer> nextList = new ArrayList<>();
							s[i] = NONE;
							List<Integer> tmp = new ArrayList<>();
							nextList.add(i);
							int p = i;
							while (true) {
								tmp.clear();
								int x = getX(p), np, dp;
								np = p + 2;
								dp = p + 1;
								if (x + 2 < N && s[dp] != NONE && s[np] == NONE) {
									tmp.add(np);
								}
								np = p - 2;
								dp = p - 1;
								if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE) {
									tmp.add(np);
								}
								np = p + N2;
								dp = p + N;
								if (p + N2 < NN && s[dp] != NONE && s[np] == NONE) {
									tmp.add(np);
								}
								np = p - N2;
								dp = p - N;
								if (p - N2 >= 0 && s[dp] != NONE && s[np] == NONE) {
									tmp.add(np);
								}
								if (tmp.isEmpty())
									break;
								np = tmp.get(rnd.nextInt(tmp.size()));
								nextList.add(np);
								s[delPos(p, np)] = NONE;
								p = np;
							}
							Collections.reverse(nextList);
							next = to(nextList);
						}
						boolean keep = true;
						Arrays.fill(index, -1);
						for (int k = 0; k < next.length; ++k)
							index[next[k]] = k;
						while (keep) {
							keep = false;
							for (int k = 1; k < next.length - 1; ++k) {
								class Inner {
									int[] func(int next[], int p, int initp) {
										Arrays.fill(prev, -1);
										prev[initp] = p;
										dist[initp] = 1;
										queue[0] = initp;
										int qs = 0, qe = 1;
										while (qs < qe) {
											int qp = queue[qs++];
											int x = getX(qp), np, dp;
											np = qp + 2;
											dp = qp + 1;
											if (x + 2 < N && s[dp] != NONE && s[np] == NONE && prev[np] == -1
													&& prev[qp] != np) {
												prev[np] = qp;
												dist[np] = dist[qp] + 1;
												if (index[np] != -1 && Math.abs(index[p] - index[np]) < dist[np]) {
													int add = dist[np] - Math.abs(index[p] - index[np]);
													int start = Math.min(index[p], index[np]);
													for (int w = start, end = Math.max(index[p], index[np]); w < end; ++w) {
														s[delPos(next[w], next[w + 1])] = 1;
													}
													int tmp[] = Arrays.copyOf(next, next.length + add);
													System.arraycopy(next, start, tmp, start + add, next.length - start);
													if (index[p] < index[np]) {
														for (int tp = np, nindex = dist[np]; nindex > 0; --nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													} else {
														for (int tp = np, nindex = 0; nindex < dist[np]; ++nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													}
													Arrays.fill(index, -1);
													for (int w = 0; w < tmp.length; ++w)
														index[tmp[w]] = w;
													return tmp;
												}
												queue[qe++] = np;
											}
											np = qp - 2;
											dp = qp - 1;
											if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE && prev[np] == -1
													&& prev[qp] != np) {
												prev[np] = qp;
												dist[np] = dist[qp] + 1;
												if (index[np] != -1 && Math.abs(index[p] - index[np]) < dist[np]) {
													int add = dist[np] - Math.abs(index[p] - index[np]);
													int start = Math.min(index[p], index[np]);
													for (int w = start, end = Math.max(index[p], index[np]); w < end; ++w) {
														s[delPos(next[w], next[w + 1])] = 1;
													}
													int tmp[] = Arrays.copyOf(next, next.length + add);
													System.arraycopy(next, start, tmp, start + add, next.length - start);
													if (index[p] < index[np]) {
														for (int tp = np, nindex = dist[np]; nindex > 0; --nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													} else {
														for (int tp = np, nindex = 0; nindex < dist[np]; ++nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													}
													Arrays.fill(index, -1);
													for (int w = 0; w < tmp.length; ++w)
														index[tmp[w]] = w;
													return tmp;
												}
												queue[qe++] = np;
											}
											np = qp + N2;
											dp = qp + N;
											if (np < NN && s[dp] != NONE && s[np] == NONE && prev[np] == -1
													&& prev[qp] != np) {
												prev[np] = qp;
												dist[np] = dist[qp] + 1;
												if (index[np] != -1 && Math.abs(index[p] - index[np]) < dist[np]) {
													int add = dist[np] - Math.abs(index[p] - index[np]);
													int start = Math.min(index[p], index[np]);
													for (int w = start, end = Math.max(index[p], index[np]); w < end; ++w) {
														s[delPos(next[w], next[w + 1])] = 1;
													}
													int tmp[] = Arrays.copyOf(next, next.length + add);
													System.arraycopy(next, start, tmp, start + add, next.length - start);
													if (index[p] < index[np]) {
														for (int tp = np, nindex = dist[np]; nindex > 0; --nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													} else {
														for (int tp = np, nindex = 0; nindex < dist[np]; ++nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													}
													Arrays.fill(index, -1);
													for (int w = 0; w < tmp.length; ++w)
														index[tmp[w]] = w;
													return tmp;
												}
												queue[qe++] = np;
											}
											np = qp - N2;
											dp = qp - N;
											if (np >= 0 && s[dp] != NONE && s[np] == NONE && prev[np] == -1
													&& prev[qp] != np) {
												prev[np] = qp;
												dist[np] = dist[qp] + 1;
												if (index[np] != -1 && Math.abs(index[p] - index[np]) < dist[np]) {
													int add = dist[np] - Math.abs(index[p] - index[np]);
													int start = Math.min(index[p], index[np]);
													for (int w = start, end = Math.max(index[p], index[np]); w < end; ++w) {
														s[delPos(next[w], next[w + 1])] = 1;
													}
													int tmp[] = Arrays.copyOf(next, next.length + add);
													System.arraycopy(next, start, tmp, start + add, next.length - start);
													if (index[p] < index[np]) {
														for (int tp = np, nindex = dist[np]; nindex > 0; --nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													} else {
														for (int tp = np, nindex = 0; nindex < dist[np]; ++nindex) {
															int prevp = prev[tp];
															s[delPos(tp, prevp)] = NONE;
															tmp[start + nindex] = tp;
															tp = prevp;
														}
													}
													Arrays.fill(index, -1);
													for (int w = 0; w < tmp.length; ++w)
														index[tmp[w]] = w;
													return tmp;
												}
												queue[qe++] = np;
											}
										}
										return null;
									}
								}
								update: {
									Inner inner = new Inner();
									int p = next[k], x = getX(p), dp, np;
									dist[p] = 0;
									dp = p + 1;
									np = p + 2;
									if (x + 2 < N && s[dp] != NONE && s[np] == NONE) {
										int tmp[] = inner.func(next, p, np);
										if (tmp != null) {
											next = tmp;
											keep = true;
											break update;
										}
									}
									dp = p - 1;
									np = p - 2;
									if (x - 2 >= 0 && s[dp] != NONE && s[np] == NONE) {
										int tmp[] = inner.func(next, p, np);
										if (tmp != null) {
											next = tmp;
											keep = true;
											break update;
										}
									}
									dp = p + N;
									np = p + N2;
									if (np < NN && s[dp] != NONE && s[np] == NONE) {
										int tmp[] = inner.func(next, p, np);
										if (tmp != null) {
											next = tmp;
											keep = true;
											break update;
										}
									}
									dp = p - N;
									np = p - N2;
									if (np >= 0 && s[dp] != NONE && s[np] == NONE) {
										int tmp[] = inner.func(next, p, np);
										if (tmp != null) {
											next = tmp;
											keep = true;
											break update;
										}
									}
								}
							}
						}
						s[next[0]] = startCell;
						if (res.score < next.length) {
							res = new State(this, next, s, next.length);
						}
					}
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

	private int[] to(List<Integer> list) {
		int res[] = new int[list.size()];
		for (int i = 0; i < res.length; ++i)
			res[i] = list.get(i);
		return res;
	}

	private static final int delPos(int p1, int p2) {
		return (p1 + p2) >> 1;
	}
}
