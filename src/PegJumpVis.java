/*
Change log
----------
2015-06-22 :
Initial release
*/

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.security.SecureRandom;

import javax.swing.JFrame;
import javax.swing.JPanel;

class TestCase {
	final Object worldLock = new Object();

	public static final int MIN_BOARD_SIZE = 20;
	public static final int MAX_BOARD_SIZE = 60;
	public static final int MIN_PEG_TYPES = 1;
	public static final int MAX_PEG_TYPES = 10;
	public static final int MIN_PEG_VALUE = 1;
	public static final int MAX_PEG_VALUE = 10;
	public static final int MIN_PEG_BLOBS = 0;
	public static final int MAX_PEG_BLOBS = 10;
	public static final double MIN_FILL_RATIO = 0.1;
	public static final double MAX_FILL_RATIO = 0.9;

	public int boardSize;

	public int pegTypeCnt;
	public int pegStartCnt;
	public int pegNowCnt;
	public int numMoves;
	public int numJumps;
	public int[] pegValue;
	public int score;

	public int blobCnt;
	public double fillRatio;

	public int[] lastX;
	public int[] lastY;
	public int[] lastPeg;
	public int lastScore;

	public char[][] board;

	public SecureRandom rnd = null;

	public TestCase(long seed) {

		try {
			rnd = SecureRandom.getInstance("SHA1PRNG");
		} catch (Exception e) {
			System.err.println("ERROR: unable to generate test case.");
			System.exit(1);
		}

		rnd.setSeed(seed);

		boardSize = rnd.nextInt(MAX_BOARD_SIZE - MIN_BOARD_SIZE + 1) + MIN_BOARD_SIZE;
		if (seed == 1)
			boardSize = 20;
		board = new char[boardSize][boardSize];
		pegTypeCnt = rnd.nextInt(MAX_PEG_TYPES - MIN_PEG_TYPES + 1) + MIN_PEG_TYPES;
		pegValue = new int[pegTypeCnt];
		for (int i = 0; i < pegTypeCnt; i++) {
			pegValue[i] = rnd.nextInt(MAX_PEG_VALUE - MIN_PEG_VALUE + 1) + MIN_PEG_VALUE;
		}
		fillRatio = rnd.nextDouble() * (MAX_FILL_RATIO - MIN_FILL_RATIO) + MIN_FILL_RATIO;
		// fill the board
		for (int y = 0; y < boardSize; y++) {
			for (int x = 0; x < boardSize; x++) {
				if (rnd.nextDouble() < fillRatio) {
					board[y][x] = '.';
				} else {
					board[y][x] = (char) ('0' + rnd.nextInt(pegTypeCnt));
				}
			}
		}
		// add the blobs
		blobCnt = rnd.nextInt(MAX_PEG_BLOBS - MIN_PEG_BLOBS + 1) + MIN_PEG_BLOBS;
		for (int b = 0; b < blobCnt; b++) {
			int pegType = rnd.nextInt(pegTypeCnt);
			int radi = rnd.nextInt(boardSize / 4) + 1;
			int centerX = rnd.nextInt(boardSize);
			int centerY = rnd.nextInt(boardSize);
			for (int x = centerX - radi; x <= centerX + radi; x++)
				if (x >= 0 && x < boardSize) {
					for (int y = centerY - radi; y <= centerY + radi; y++)
						if (y >= 0 && y < boardSize) {
							int r = (x - centerX) * (x - centerX) + (y - centerY) * (y - centerY);
							if (r <= radi * radi)
								board[y][x] = (char) ('0' + pegType);
						}
				}
		}
		pegStartCnt = 0;
		for (int y = 0; y < boardSize; y++) {
			for (int x = 0; x < boardSize; x++) {
				if (board[y][x] != '.')
					pegStartCnt++;
			}
		}
		pegNowCnt = pegStartCnt;
		numJumps = 0;
		score = 0;
		lastScore = -1;
	}

	public boolean doMove(int x, int y, String mv) throws Exception {
		synchronized (worldLock) {
			char peg = board[y][x];
			if (peg == '.') {
				System.err.println("ERROR: Trying to move an empty cell.");
				return false;
			}
			lastX = new int[mv.length() + 1];
			lastY = new int[mv.length() + 1];
			lastPeg = new int[mv.length() + 1];
			lastX[0] = x;
			lastY[0] = y;
			int sum = 0;
			for (int i = 0; i < mv.length(); i++) {
				int dx = 0;
				int dy = 0;
				char chMove = mv.charAt(i);
				if (chMove == 'L')
					dx = -1;
				else if (chMove == 'R')
					dx = 1;
				else if (chMove == 'U')
					dy = -1;
				else if (chMove == 'D')
					dy = 1;
				if (dx == 0 && dy == 0) {
					System.err.println("ERROR: Invalid move character [" + chMove + "].");
					return false;
				}
				int nx = x + dx;
				int ny = y + dy;
				if (nx >= 0 && nx < boardSize && ny >= 0 && ny < boardSize) {
					if (board[ny][nx] == '.') {
						System.err.println("ERROR: Can not jump over empty space.");
						System.err.println(mv);
						return false;
					} else {
						int nx2 = nx + dx;
						int ny2 = ny + dy;
						if (nx2 >= 0 && nx2 < boardSize && ny2 >= 0 && ny2 < boardSize) {
							if (board[ny2][nx2] != '.') {
								System.err.println("ERROR: Trying to jump onto another peg.");
								return false;
							}
							lastPeg[i + 1] = (board[ny][nx] - '0');
							board[y][x] = '.';
							sum += pegValue[(board[ny][nx] - '0')];
							board[ny][nx] = '.';
							board[ny2][nx2] = peg;
							x = nx2;
							y = ny2;
							lastX[i + 1] = x;
							lastY[i + 1] = y;
							pegNowCnt--;
						} else {
							System.err.println("ERROR: Trying to jump outside the board.");
							return false;
						}
					}
				} else {
					System.err.println("ERROR: Trying to move outside the board.");
					return false;
				}
			}
			lastScore = mv.length() * sum;
			score += lastScore;
			numJumps++;
		}
		return true;
	}
}

class Drawer extends JFrame {
	public static final int EXTRA_WIDTH = 200;
	public static final int EXTRA_HEIGHT = 50;

	public TestCase tc;
	public DrawerPanel panel;

	public int cellSize, boardSize;
	public int width, height;

	public boolean pauseMode = false;
	public boolean debugMode = false;

	class DrawerKeyListener extends KeyAdapter {
		public void keyPressed(KeyEvent e) {
			synchronized (keyMutex) {
				if (e.getKeyChar() == ' ') {
					pauseMode = !pauseMode;
				}
				if (e.getKeyChar() == 'd') {
					debugMode = !debugMode;
				}
				keyPressed = true;
				keyMutex.notifyAll();
			}
		}
	}

	class DrawerPanel extends JPanel {
		public void paint(Graphics g) {
			synchronized (tc.worldLock) {
				g.setColor(new Color(32, 32, 32));
				g.fillRect(15, 15, cellSize * boardSize + 1, cellSize * boardSize + 1);
				g.setColor(Color.BLACK);
				for (int i = 0; i <= boardSize; i++) {
					g.drawLine(15 + i * cellSize, 15, 15 + i * cellSize, 15 + cellSize * boardSize);
					g.drawLine(15, 15 + i * cellSize, 15 + cellSize * boardSize, 15 + i * cellSize);
				}
				for (int i = 0; i < boardSize; i++) {
					for (int j = 0; j < boardSize; j++) {
						if (tc.board[i][j] >= '0' && tc.board[i][j] <= '9') {
							// peg color
							float hue = (float) (tc.board[i][j] - '0') / tc.pegTypeCnt;
							g.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
							g.fillOval(15 + j * cellSize + 1, 15 + i * cellSize + 1, cellSize - 1, cellSize - 1);
						}
					}
				}

				if (tc.lastX != null && tc.lastX.length > 0) {
					for (int i = 0; i < tc.lastX.length - 1; i++) {

						float hue = (float) (tc.lastPeg[i + 1]) / tc.pegTypeCnt;
						g.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
						g.drawOval(15 + ((tc.lastX[i + 1] + tc.lastX[i]) / 2) * cellSize + 1, 15
								+ ((tc.lastY[i + 1] + tc.lastY[i]) / 2) * cellSize + 1, cellSize - 1, cellSize - 1);

						g.setColor(Color.WHITE);
						float f = 0.5f + 0.5f * (i + 1) / (tc.lastX.length - 1);
						g.setColor(new Color(f, f, f));

						g.drawLine(15 + tc.lastX[i] * cellSize + cellSize / 2, 15 + tc.lastY[i] * cellSize + cellSize
								/ 2, 15 + tc.lastX[i + 1] * cellSize + cellSize / 2, 15 + tc.lastY[i + 1] * cellSize
								+ cellSize / 2);
						g.fillOval(15 + tc.lastX[i] * cellSize + cellSize / 4, 15 + tc.lastY[i] * cellSize + cellSize
								/ 4, cellSize / 2, cellSize / 2);
					}
				}
				int horPos = 40 + boardSize * cellSize;
				g.setColor(Color.BLACK);
				g.setFont(new Font("Arial", Font.BOLD, 12));
				Graphics2D g2 = (Graphics2D) g;
				g2.drawString("Board size = " + boardSize, horPos, 30);
				g2.drawString("Pegs at start = " + tc.pegStartCnt, horPos, 50);
				g2.drawString("Pegs on board = " + tc.pegNowCnt, horPos, 70);
				g2.drawString("Jumps = " + tc.numJumps, horPos, 90);
				if (tc.lastScore > 0)
					g2.drawString("Sequence score = " + tc.lastScore, horPos, 110);
				g2.drawString("Score = " + tc.score, horPos, 140);
				for (int i = 0; i < tc.pegTypeCnt; i++) {
					float hue = (float) (i) / tc.pegTypeCnt;
					g.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
					g.fillOval(horPos, 170 + i * cellSize + 1, cellSize - 1, cellSize - 1);
					g.setColor(Color.BLACK);
					g2.drawString(Integer.toString(tc.pegValue[i]), horPos + cellSize * 2, 180 + i * cellSize + 1);
				}

			}
		}
	}

	class DrawerWindowListener extends WindowAdapter {
		public void windowClosing(WindowEvent event) {
			System.exit(0);
		}
	}

	final Object keyMutex = new Object();
	boolean keyPressed;

	public void processPause() {
		synchronized (keyMutex) {
			if (!pauseMode) {
				return;
			}
			keyPressed = false;
			while (!keyPressed) {
				try {
					keyMutex.wait();
				} catch (InterruptedException e) {
					// do nothing
				}
			}
		}
	}

	public Drawer(TestCase tc_, int cellSize) {
		super();

		panel = new DrawerPanel();
		getContentPane().add(panel);

		addWindowListener(new DrawerWindowListener());

		this.tc = tc_;

		boardSize = tc.boardSize;
		this.cellSize = cellSize;
		width = cellSize * boardSize + EXTRA_WIDTH;
		height = cellSize * boardSize + EXTRA_HEIGHT;

		addKeyListener(new DrawerKeyListener());

		setSize(width, height);
		setTitle("Visualizer tool for problem PegJump");

		setResizable(false);
		setVisible(true);
	}
}

public class PegJumpVis {
	public static boolean vis = true;
	public static boolean debug = false;
	public static int cellSize = 12;
	public static int delay = 100;
	public static boolean startPaused = false;

	public int runTest(long seed) {

		TestCase tc = new TestCase(seed);

		Drawer drawer = null;
		if (vis) {
			drawer = new Drawer(tc, cellSize);
			drawer.debugMode = debug;
			if (startPaused) {
				drawer.pauseMode = true;
			}
		}

		try {
			String board[] = new String[tc.board.length];
			for (int i = 0; i < tc.board.length; ++i) {
				board[i] = new String(tc.board[i]);
			}
			long start = System.currentTimeMillis();
			String[] res = new PegJumping().getMoves(tc.pegValue, board);
			System.out.println("Time = " + (System.currentTimeMillis() - start));
			int numMoves = res.length;
			if (numMoves > tc.boardSize * tc.boardSize * 10) {
				System.err.println("ERROR: Return array from getMoves too large.");
				return -1;
			}
			for (int i = 0; i < numMoves; i++) {
				String smove = res[i];
				String[] s = smove.split(" ");
				if (s.length != 3) {
					System.err.println("ERROR: The move command with index " + i
							+ " does not contain 3 space seperated values. Value is [" + smove + "]");
					return -1;
				}
				int y = Integer.parseInt(s[0]);
				int x = Integer.parseInt(s[1]);
				if (x < 0 || x >= tc.boardSize || y < 0 || y >= tc.boardSize) {
					System.err.println("ERROR: (" + y + "," + x + ") outside of bounds.");
					return -1;
				}
				if (tc.board[y][x] == '.') {
					System.err.println("ERROR: (" + y + "," + x + ") does not contain a peg.");
					return -1;
				}
				// perform moves
				if (!tc.doMove(x, y, s[2])) {
					// error occured
					return -1;
				}
				if (vis) {
					drawer.processPause();
					drawer.repaint();
					try {
						Thread.sleep(delay);
					} catch (Exception e) {
						// do nothing
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("ERROR: Unable to process the move commands from your solution.");
			return -1;
		}
		return tc.score;
	}

	public static void main(String[] args) {
		for (int i = 0; i < args.length; i++)
			if (args[i].equals("-novis")) {
				vis = false;
			} else if (args[i].equals("-debug")) {
				debug = true;
			} else if (args[i].equals("-sz")) {
				cellSize = Integer.parseInt(args[++i]);
			} else if (args[i].equals("-delay")) {
				delay = Integer.parseInt(args[++i]);
			} else if (args[i].equals("-pause")) {
				startPaused = true;
			} else {
				System.out.println("WARNING: unknown argument " + args[i] + ".");
			}

		vis = false;
		PegJumpVis vis = new PegJumpVis();
		try {
			for (long seed = 1; seed <= 1000; ++seed) {
				int score = vis.runTest(seed);
				System.out.println("Seed = " + seed + "   Score = " + score);
			}
		} catch (RuntimeException e) {
			System.err.println("ERROR: Unexpected error while running your test case.");
			e.printStackTrace();
		}
	}
}

class ErrorStreamRedirector extends Thread {
	public BufferedReader reader;

	public ErrorStreamRedirector(InputStream is) {
		reader = new BufferedReader(new InputStreamReader(is));
	}

	public void run() {
		while (true) {
			String s;
			try {
				s = reader.readLine();
			} catch (Exception e) {
				//e.printStackTrace();
				return;
			}
			if (s == null) {
				break;
			}
			System.err.println(s);
		}
	}
}