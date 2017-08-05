#pragma GCC optimize "O3,omit-frame-pointer,inline,fast-math"

//collision algorithm based on http://files.magusgeek.com/csb/csb_en.html

#include <iostream>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <memory>
#include <chrono>

#define PI 3.14159265358979323846
//default max time allowed for IA computing (in ms)
#define DEFAULT_MAX_TIME 145
//number of moves per pod in a solution
#define LENGTH 6
//max distance on map, used for normalization
#define MAX_DIST 1900.0
//shield activation probability
#define ASHIELD 0.1
//probability for max and min thrust
#define MAXTHRUSTP 0.1
#define MINTHRUSTP 0.1


//exponential function approximation
inline
float exp(float x) {
	  x = 1.0 + x / 16.0;
	  x *= x; x *= x; x *= x; x *= x;
	  return x;
}

//singleton class to handle random numbers
class Random {
private:
    std::mt19937 engine ;
    std::uniform_real_distribution<float> fdistribution ;
    std::uniform_int_distribution<int> angle_distribution ;
    std::uniform_int_distribution<int> thrust_distribution ;
    std::uniform_int_distribution<int> length_distribution ;

    //copies are not allowed
    Random(Random const&) {}
	void operator=(Random const&) {}

	//constructor is private
	//getInstance is the way to get a random object (unique)
	Random() {
	    static std::random_device rd; // obtain a random number from hardware
	    engine = std::mt19937(rd()); // seed the generator
	    fdistribution = std::uniform_real_distribution<float>(0.0,1.0) ;
	    angle_distribution = std::uniform_int_distribution<int>(-18,18) ;
	    thrust_distribution = std::uniform_int_distribution<int>(0, 200) ;
	    length_distribution = std::uniform_int_distribution<int>(0, 2 * LENGTH - 1) ;
	}

public:
	~Random(){}

	static Random& getInstance() {
		static Random instance ;
		return instance ;
	}
	float nextFloat() {
		return fdistribution(engine) ;
	}
	int nextRAngle() {
		return angle_distribution(engine);
	}
	int nextRThrust() {
		return thrust_distribution(engine);
	}
	int nextRLength() {
		return length_distribution(engine);
	}
};

class Move {
public:
	Move(float angle = 0., int thrust=0){
		this->angle = angle ;
		this->thrust = thrust ;
	}
	Move(Move const &m) {
		this->angle = m.angle ;
		this->thrust = m.thrust ;
	}
	float getAngle() const {
		return angle;
	}

	void setAngle(float angle) {
		this->angle = angle;
	}

	int getThrust() const {
		return thrust;
	}

	void setThrust(int thrust) {
		this->thrust = thrust;
	}

	void mutate() {
		Random &r = Random::getInstance() ;
		if(r.nextFloat() < ASHIELD) {
			this->thrust = this->thrust == -1 ? r.nextRThrust() : -1 ;
		}
		else if(r.nextFloat() < ASHIELD + MAXTHRUSTP)
		{
			this->thrust = 200;
		}
		else if(r.nextFloat() < ASHIELD + MAXTHRUSTP + MINTHRUSTP) {
			this->thrust = 0 ;
		}
		else {
			this->thrust = r.nextRThrust() ;
		}
		this->angle = r.nextRAngle() ;
	}

	static Move randomMove() {
		Random &r = Random::getInstance() ;
		float thrust ;
		if(r.nextFloat() < ASHIELD) thrust = -1 ;
		else thrust = r.nextRThrust() ;

		return Move(r.nextRAngle(), thrust) ;
	}

private :
	float angle; // Between -18.0 and +18.0
    int thrust; // Between -1 and 200
};

std::ostream& operator <<(std::ostream& Stream, const Move& move)
{
    Stream << "Angle : " << move.getAngle() << " Thrust : " << move.getThrust() ;
    return Stream; // N'oubliez pas de renvoyer le flux, afin de pouvoir chaîner les appels
}

class Solution {
private:
	float score ;

public:
	//LENGTH first elements : first pod
	//LENGTH last elements : second pod
	Move solution[2 * LENGTH] ;

	Solution(bool buildRandom = false) {
		score = std::numeric_limits<float>::lowest();
		if(buildRandom) {
			for(int i = 0 ; i < 2* LENGTH ; i++) {
				this->solution[i] = Move::randomMove() ;
			}
		}
		else {
			for(int i = 0 ; i < 2* LENGTH ; i++) {
				this->solution[i] = Move() ;
			}
		}
	}

	Solution(const Solution& s) {
		this->score = s.score ;
		std::copy(std::begin(s.solution), std::end(s.solution), std::begin(this->solution)) ;
	}

	void shiftLeft() {
		std::copy(std::begin(this->solution)+1, std::end(this->solution),std::begin(this->solution)) ;
		this->solution[LENGTH - 1] = Move(0,0) ;
	}

	static Solution* randomSolution() {
		Solution* s = new Solution() ;
		for(int i = 0 ; i < 2* LENGTH ; i++) {
			s->solution[i] = Move::randomMove() ;
		}
		return s ;
	}

	void mutate(){
		Random& r = Random::getInstance() ;
		solution[r.nextRLength()].mutate() ;
	}

	float getScore() const {
		return score;
	}

	void setScore(float score) {
		this->score = score;
	}
};

class Point {
public:
	float x ;
	float y ;

	Point() {
		x = 0. ;
		y = 0. ;
	}
	Point(float x, float y) {
		this->x = x ;
		this->y = y ;
	}
	Point(Point const &p) {
		this->x = p.x ;
		this->y = p.y ;
	}
	float distance(float x, float y) {
		return sqrt((this->x - x) * (this->x - x) + (this->y - y) * (this->y - y)) ;
	}
	float distance(const Point &p) {
		return this->distance(p.x, p.y) ;
	}
	float distanceSq(float x, float y) {
		return (this->x - x) * (this->x - x) + (this->y - y) * (this->y - y) ;
	}
	float distanceSq(const Point &p) {
		return this->distanceSq(p.x, p.y) ;
	}

	float getX() const {
		return x;
	}

	void setX(float x) {
		this->x = x;
	}

	float getY() const {
		return y;
	}

	void setY(float y) {
		this->y = y;
	}

	//It allows us to find the closest point on a line (described by two points) to another point.
	Point closest(const Point &a, const Point &b) {
	    float da = b.y - a.y;
	    float db = a.x - b.x;
	    float c1 = da*a.x + db*a.y;
	    float c2 = -db*this->x + da*this->y;
	    float det = da*da + db*db;
	    float cx = 0;
	    float cy = 0;

	    if (det != 0) {
	        cx = (da*c1 - db*c2) / det;
	        cy = (da*c2 + db*c1) / det;
	    } else {
	        // The point is already on the line
	        cx = this->x;
	        cy = this->y;
	    }

	    return Point(cx, cy);
	}

protected:

};

//forward declarations
class Unit ;
class Pod ;

class Collision {
public:
	Unit* a ;
	Unit* b;
	float t ;

	Collision(Unit* a, Unit* b, float t) :
		a(a), b(b), t(t) {

	}

private:

};

class Unit:public Point {
public:
	Unit(float x, float y):
		Point(x, y)
	{
		r = 0 ;
		vx = 0 ;
		vy = 0 ;
	}

	virtual ~Unit() {
	}

	Unit(float x, float y, float vx, float vy):
		Point(x, y)
	{
		r = 0 ;
		this->vx = vx ;
		this->vy = vy ;
	}

	float getR() const {
		return r;
	}

	float getVx() const {
		return vx;
	}

	void setVx(float vx) {
		this->vx = vx;
	}

	float getVy() const {
		return vy;
	}

	void setVy(float vy) {
		this->vy = vy;
	}

	virtual void bounce(Unit& unit) {
		throw std::logic_error("Appel à la fonction virtuelle Bounce de la classe Unit") ;
	}

	Collision* collision(Unit* u) {
		// Square of the distance
	    float dist = this->distanceSq(*u);

	    // Sum of the radii squared
	    float sr = (this->r + u->r)*(this->r + u->r);

	    // We take everything squared to avoid calling sqrt uselessly. It is better for performances

	    if (dist < sr) {
	        // Objects are already touching each other. We have an immediate collision.
	        return new Collision(this, u, 0.0);
	    }

	    // Optimisation. Objects with the same speed will never collide
	    if (this->vx == u->vx && this->vy == u->vy) {
	        return NULL;
	    }

	    // We place ourselves in the reference frame of u. u is therefore stationary and is at (0,0)
	    float x = this->x - u->x;
	    float y = this->y - u->y;
	    Point myp(x, y);
	    float vx = this->vx - u->vx;
	    float vy = this->vy - u->vy;
	    Point up(0., 0.) ;

	    // We look for the closest point to u (which is in (0,0)) on the line described by our speed vector
	    Point p = up.closest(myp, Point(x + vx, y + vy));

	    // Square of the distance between u and the closest point to u on the line described by our speed vector
	    float pdist = up.distanceSq(p);

	    // Square of the distance between us and that point
	    float mypdist = myp.distanceSq(p);

	    // If the distance between u and this line is less than the sum of the radii, there might be a collision
	    if (pdist < sr) {
	     // Our speed on the line
	        float length = sqrt(vx*vx + vy*vy);

	        // We move along the line to find the point of impact
	        float backdist = sqrt(sr - pdist);
	        p.x = p.x - backdist * (vx / length);
	        p.y = p.y - backdist * (vy / length);

	        // If the point is now further away it means we are not going the right way, therefore the collision won't happen
	        if (myp.distanceSq(p) > mypdist) {
	            return NULL;
	        }

	        pdist = p.distance(myp);

	        // The point of impact is further than what we can travel in one turn
	        if (pdist > length) {
	            return NULL;
	        }

	        // Time needed to reach the impact point
	        float t = pdist / length;

	        return new Collision(this, u, t);
	    }

	    return NULL;
	}


protected:
	float r ;
	float vx ;
	float vy ;
} ;



class Pod: public Unit {
public:
	Pod(float x = -1., float y=-1., float vx=0., float vy=0., float angle = 0., int nextCheckpointId = -1):
		Unit(x,y,vx,vy)
	{
		this->r = 400 ;
		this->angle = angle ;
		this->nextCheckpointId = nextCheckpointId ;
		this->timeout = 100 ;
		this->shield = 0 ;
		this->checked = 0 ;
	}

	float getAngle() const {
		return angle;
	}

	void setAngle(float angle) {
		this->angle = angle;
	}

	int getNextCheckpointId() const {
		return nextCheckpointId;
	}

	void setNextCheckpointId(int nextCheckpointId) {
		this->nextCheckpointId = nextCheckpointId;
	}

	bool isShieldActive() {
		return this->shield == 3 ;
	}

	int getTimeout() {
		return this->timeout ;
	}

	void checkedCP() {
		this->timeout = 100 ;
		this->checked ++ ;
		this->nextCheckpointId ++ ;
	}

	//compute the angle that we should have to target a point
	float computeAngle(const Point& p) {
	    float d = this->distance(p);
	    float dx = (p.x - this->x) / d;
	    float dy = (p.y - this->y) / d;

	    // Simple trigonometry. We multiply by 180.0 / PI to convert radiants to degrees.
	    float a = acos(dx) * 180.0 / PI;

	    // If the point I want is below me, I have to shift the angle for it to be correct
	    if (dy < 0) {
	        a = 360.0 - a;
	    }
	    return a;
	}

	//return the difference between the target angle and the current angle
	float diffAngle(const Point& p) {
	    float a = this->computeAngle(p);

	    // To know whether we should turn clockwise or not we look at the two ways and keep the smallest
	    // The ternary operators replace the use of a modulo operator which would be slower
	    float right = this->angle <= a ? a - this->angle : 360.0 - this->angle + a;
	    float left = this->angle >= a ? this->angle - a : this->angle + 360.0 - a;

	    if (right < left) {
	        return right;
	    } else {
	        // We return a negative angle if we must rotate to left
	        return -left;
	    }
	}

	float diffSpeedAngle(const Point& p) {
		if(this->vx == 0 && this->vy == 0) return 0. ;
		//target angle
		float a = this->computeAngle(p) ;
		//angle given by speed vector
		float sa = this->computeAngle(Point(x + this->vx, y + this->vy)) ;

	    // To know whether we should turn clockwise or not we look at the two ways and keep the smallest
	    // The ternary operators replace the use of a modulo operator which would be slower
	    float right = sa <= a ? a - sa : 360.0 - sa + a;
	    float left = sa >= a ? sa - a : sa + 360.0 - a;

	    if (right < left) {
	        return right;
	    } else {
	        // We return a negative angle if we must rotate to left
	        return -left;
	    }
	}

	void rotate(const Point& p) {
	    float a = this->diffAngle(p);

	    // Can't turn by more than 18° in one turn
	    if (a > 18.0) {
	        a = 18.0;
	    } else if (a < -18.0) {
	        a = -18.0;
	    }

	    this->angle += a;

	    // The % operator is slow. If we can avoid it, it's better.
	    if (this->angle >= 360.0) {
	        this->angle = this->angle - 360.0;
	    } else if (this->angle < 0.0) {
	        this->angle += 360.0;
	    }
	}

	void rotate(float angle) {
		if(angle > 18.0) angle = 18.0 ;
		else if(angle < -18.0) angle = -18.0 ;

		this->angle += angle;

		// The % operator is slow. If we can avoid it, it's better.
		if (this->angle >= 360.0) {
		    this->angle = this->angle - 360.0;
		} else if (this->angle < 0.0) {
		    this->angle += 360.0;
		}
	}

	void boost(int thrust) {
	    if(thrust == -1) {
	    	this->shield = 3 ;
	    	return;
	    }
		//a pod which has activated its shield cannot accelerate for 3 turns
	    if (this->shield > 0) {
	        return;
	    }
	    // Conversion of the angle to radiants
	    float ra = this->angle * PI / 180.0;

	    // Trigonometry
	    this->vx += cos(ra) * thrust;
	    this->vy += sin(ra) * thrust;
	}

	void move(float t) {
	    this->x += this->vx * t;
	    this->y += this->vy * t;
	}

	void endTurn(int checkpointCount) {
	    this->x = round(this->x);
	    this->y = round(this->y);
	    this->vx = (int)(this->vx * 0.85);
	    this->vy = (int)(this->vy * 0.85);

	    // Don't forget that the timeout goes down by 1 each turn. It is reset to 100 when you pass a checkpoint
	    this->timeout -= 1;
	    if(this->shield > 0) this->shield -- ;
		if(this->nextCheckpointId == checkpointCount) this->nextCheckpointId = 0 ;
	}

	void bounce(Unit& unit) {
		Pod& u = (Pod&)unit ;

		// If a pod has its shield active its mass is 10 otherwise it's 1
		float m1 = this->isShieldActive() ? 10 : 1;
		float m2 = u.isShieldActive() ? 10 : 1;
		float mcoeff = (m1 + m2) / (m1 * m2);

		float nx = this->x - u.x;
		float ny = this->y - u.y;

		// Square of the distance between the 2 pods. This value could be hardcoded because it is always 800²
		//TODO
		float nxnysquare = nx*nx + ny*ny;

		float dvx = this->vx - u.vx;
		float dvy = this->vy - u.vy;

		// fx and fy are the components of the impact vector. product is just there for optimisation purposes
		float product = nx*dvx + ny*dvy;
		float fx = (nx * product) / (nxnysquare * mcoeff);
		float fy = (ny * product) / (nxnysquare * mcoeff);

		// We apply the impact vector once
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u.vx += fx / m2;
		u.vy += fy / m2;

		// If the norm of the impact vector is less than 120, we normalize it to 120
		float impulse = sqrt(fx*fx + fy*fy);
		if (impulse < 120.0) {
			fx = fx * 120.0 / impulse;
			fy = fy * 120.0 / impulse;
		}

		// We apply the impact vector a second time
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u.vx += fx / m2;
		u.vy += fy / m2;
	}

	void Output(Move& move) {
	    std::string s ;
		float a = angle + move.getAngle();

	    if (a >= 360.0) {
	        a = a - 360.0;
	    } else if (a < 0.0) {
	        a += 360.0;
	    }

	    // Look for a point corresponding to the angle we want
	    // Multiply by 10000.0 to limit rounding errors
	    a = a * PI / 180.0;
	    float px = this->x + cos(a) * 10000.0;
	    float py = this->y + sin(a) * 10000.0;

	    if (move.getThrust() == -1) {
	        std::cout << round(px) << " " << round(py) << " SHIELD" << std::endl;
	    } else {
	        std::cout << round(px) << " " << round(py) <<" " << move.getThrust() << std::endl;
	    }
	}

	int getChecked() const {
		return checked;
	}

	float score(Point* dest) {
		return this->getChecked() * 100 - this->distance(*dest) / MAX_DIST ;
	}

protected:
	float angle ;
	int nextCheckpointId ;
	int timeout ;
	int shield ;
	int checked ;
};

class Checkpoint: public Unit {
public:
	Checkpoint(float x = -1., float y = -1.):
		Unit(x, y)
	{
		this->r = 200 ; //600 - 400
	}

	void bounce(Unit& unit) {
		((Pod&)unit).checkedCP() ;
	}
};

class Game {
private:
	static int laps ;
	static int checkpointCount ;

public:
	int turn ;
	static Checkpoint checkpoints[8] ;
	Pod pods[4] ;

	Game(int laps,int checkpointCount) {
	    turn = 0 ;
        Game::laps = laps ;
        Game::checkpointCount = checkpointCount ;

		for (int i = 0; i < checkpointCount; i++) {
	    	int checkpointX;
	        int checkpointY;
	        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
	        checkpoints[i] = Checkpoint(checkpointX, checkpointY) ;
	    }
	}

	Game(const Game& g) {
		turn = g.turn ;

		for(int i = 0 ; i < 4 ; i++) {
			pods[i] = Pod(g.pods[i]) ;
		}
	}

	void updatePods(int i, float x, float y, float vx, float vy, float angle, int nextCheckPointId) {
		if(turn == 1) {
			//on first turn, pad angle is undefined (-1) and calculated based on first direction
			//to avoid having a different behavior on first turn, we set the angle to next CP
            pods[i] = Pod(x,y,vx,vy,angle,nextCheckPointId) ;
			float startAngle = pods[i].computeAngle(checkpoints[pods[i].getNextCheckpointId()]) ;
			pods[i].setAngle(startAngle) ;
		}
		else {
			Pod& pod = pods[i] ;
			pod.setX(x);
			pod.setY(y);
			pod.setVx(vx);
			pod.setVy(vy);
			pod.setAngle(angle);
			if(nextCheckPointId != pod.getNextCheckpointId()) pod.checkedCP() ;
			pod.setNextCheckpointId(nextCheckPointId);
		}
	}

	void readGame() {
		turn ++;

        for (int i = 0; i < 4; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            std::cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; std::cin.ignore();
            updatePods(i, x, y, vx, vy, angle, nextCheckPointId);
        }
	}

	static const int getCheckpointCount() {
		return checkpointCount;
	}

	static const int getLaps() {
		return laps;
	}

private:
	void simulateMovement() {
	    // This tracks the time during the turn. The goal is to reach 1.0
	    float t = 0.0;
	    Collision* previousCollision = NULL ;

	    while (t < 1.0) {
	        Collision* firstCollision = NULL ;

	        // We look for all the collisions that are going to occur during the turn
	        for (int i = 0; i < 4 ; i++) {
	            // Collision with another pod?
	            for (int j = i + 1 ; j < 4 ; j++) {
	                Collision* col = pods[i].collision(&pods[j]);

	                //already played collision ?
	                //necessary t avoid infinite loops due to rounding errors
	                //TODO : fix if t = 0
	                if(col != NULL && col->t == 0) continue ;

	                if(previousCollision != NULL && col != NULL && col->t == 0 &&
	                		col->a == previousCollision->a && col->b == previousCollision->b) {
	                	delete col ;
	                	continue ;
	                }
	                // If the collision occurs earlier than the one we currently have we keep it
	                else if (col != NULL && col->t + t < 1.0 && (firstCollision == NULL || col->t < firstCollision->t)) {
	                    firstCollision = col;
	                }
	                else if(col != NULL) {
	                	delete col ;
	                }
	            }

	            // Collision with another checkpoint?
	            // It is unnecessary to check all checkpoints here. We only test the pod's next checkpoint.
	            // We could look for the collisions of the pod with all the checkpoints, but if such a collision happens it wouldn't impact the game in any way
	            Collision* col = pods[i].collision(checkpoints + pods[i].getNextCheckpointId());

	            // If the collision happens earlier than the current one we keep it
	            if (col != NULL && col->t + t < 1.0 && (firstCollision == NULL || col->t < firstCollision->t)) {
	                firstCollision = col;
	            }
	            else if(col !=NULL){
	            	delete col ;
	            }
	        }

	        if (firstCollision == NULL) {
	            // No collision, we can move the pods until the end of the turn
	            for (int i = 0; i < 4 ; i++) {
	                pods[i].move(1.0 - t);
	            }

	            // End of the turn
	            t = 1.0;
	        } else {
	            // Move the pods to reach the time `t` of the collision
	            for (int i = 0; i < 4 ; ++i) {
	                pods[i].move(firstCollision->t);
	            }

	            // Play out the collision
	            firstCollision->b->bounce(*firstCollision->a) ;
	            t += firstCollision->t;

	            delete previousCollision ;
	            previousCollision = firstCollision ;
	        }
	    }

	    delete previousCollision ;

	    for (int i = 0; i < 4 ; ++i) {
	        pods[i].endTurn(Game::getCheckpointCount());
	    }

	    turn ++ ;
	}

public:
	/*void simulateNextTurn(const std::vector<Move>& moves) {
		for(int i = 0; i < static_cast<int>(moves.size()) ; i++) {
			this->pods[i].rotate(moves[i].getAngle());
			this->pods[i].boost(moves[i].getThrust()) ;
		}
		this->simulateMovement();
	}*/

	void simulateNextTurn(const Move& move0, const Move& move1, const Move& move2, const Move& move3) {
		this->pods[0].rotate(move0.getAngle());
		this->pods[0].boost(move0.getThrust()) ;

		this->pods[1].rotate(move1.getAngle());
		this->pods[1].boost(move1.getThrust()) ;

		this->pods[2].rotate(move2.getAngle());
		this->pods[2].boost(move2.getThrust()) ;

		this->pods[3].rotate(move3.getAngle());
		this->pods[3].boost(move3.getThrust()) ;

		this->simulateMovement();
	}

	void simulateSolutions(Solution& s, Solution& os) {
		for(int i = 0 ; i < LENGTH ; i ++) {
			simulateNextTurn(s.solution[i], s.solution[i+LENGTH], os.solution[i], os.solution[i+LENGTH]) ;
		}
	}

	float evalGame(bool isOpp=false) {

		int fpi = 0, spi = 1, ofpi = 2, ospi =3 ;

		if(isOpp){
			fpi = 2;
			spi = 3 ;
			ofpi = 0 ;
			ospi = 1 ;
		}

		float score0 = pods[fpi].score(checkpoints + pods[fpi].getNextCheckpointId()) ;
		float score1 = pods[spi].score(checkpoints + pods[spi].getNextCheckpointId()) ;
		float score2 = pods[ofpi].score(checkpoints + pods[ofpi].getNextCheckpointId()) ;
		float score3 = pods[ospi].score(checkpoints + pods[ospi].getNextCheckpointId()) ;

		return score0 + score1 - std::max(score2, score3) ;
	}
};

//static variables initialization
int Game::laps = -1 ;
int Game::checkpointCount = -1 ;
Checkpoint Game::checkpoints[8] ;

class IA {
public:
	Game* game ;
	std::chrono::time_point<std::chrono::steady_clock> begin ;

	//virtual Solution* computeSolution() = 0 ;
	virtual ~IA() { } ;

	IA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp = false) :
		game(game), MAX_TIME(max_time), isOpp(isOpp) {
		begin = std::chrono::steady_clock::now();
	}

	void resetTimer() {
		begin = std::chrono::steady_clock::now();
	}

	float getElapsedTime() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count();
	}

	bool timeOut() {
		return getElapsedTime() >= MAX_TIME ;
	}

	/*template<class IA>
	static Solution* generateSolution(const Game& game) {
		Solution* s = new Solution() ;

		Game simulation = Game(game) ;
		IA ia = IA(&simulation, false) ;

		for(int i = 0 ; i < LENGTH ; i++) {
			std::vector<Move> m = ia.computeMoves() ;
			s->solution[i] = m[0] ;
			s->solution[i+LENGTH] = m[1] ;
			simulation.simulateNextTurn(m) ;
		}
		return s;
	}*/

	template<class IA>
	float scoreSolution(Solution& s) {
		Game simulation = Game(*game) ;

		IA ia = IA(&simulation, -1, isOpp) ;
		std::unique_ptr<Move[]> om = nullptr;

		//si on ne joue pas l'adversaire
		if(!isOpp) {
			for(int i = 0 ; i < LENGTH ; i++) {
				om.reset(ia.computeMoves()) ;
				simulation.simulateNextTurn(s.solution[i], s.solution[i+LENGTH], om[0], om[1]) ;
			}
		}
		//si on joue l'adversaire : l'IA joue notre rôle, la solution est appliquée à l'adv
		else {
			for(int i = 0 ; i < LENGTH ; i++) {
				om.reset(ia.computeMoves()) ;
				simulation.simulateNextTurn(om[0], om[1],s.solution[i], s.solution[i+LENGTH]) ;
			}
		}
		s.setScore(simulation.evalGame(isOpp)) ;
		return s.getScore() ;
	}

	float scoreSolution(Solution* s, Solution* os) {
		Game simulation = Game(*game) ;

		for(int i = 0 ; i < LENGTH ; i++) {
			simulation.simulateNextTurn(s->solution[i], s->solution[i+LENGTH], os->solution[i], os->solution[i+LENGTH]) ;
		}
		s->setScore(simulation.evalGame(isOpp)) ;
		return s->getScore() ;
	}

protected :
	int MAX_TIME ;
	bool isOpp ;
};

class SimpleIA : public IA {
public:
	SimpleIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp=false) : IA(game,max_time,isOpp) {
		if(!isOpp) {
			myPods = game->pods ;
			oppPods = game->pods + 2 ;
		}
		else {
			myPods = game->pods +2 ;
			oppPods = game -> pods ;
		}
	}

	~SimpleIA() {

	}
protected:
	//advance to next checkpoint
	Move computeAMove(Pod* pod) {
		Unit* target = game->checkpoints + pod->getNextCheckpointId() ;

		float speed = sqrt(pod->getVx() * pod->getVx() + pod->getVy() * pod->getVy()) ;
		float thrust = 200 ;
		float n = pod->distance(*target) / speed ;

		//if we are close to the checkpoint, we try to go to the next one
		//this allows to turn before reaching checkpoint
		if(n <= 5) {
			//new target : next checkpoint
			Unit* newTarget = game->checkpoints + (pod->getNextCheckpointId() + 1) % game->getCheckpointCount() ;
			Collision* col = NULL ;
			for(int tthrust = 200 ; tthrust >= 0 && col == NULL; tthrust-=10) {
				Pod clone = Pod(*pod) ;
				for(int j = 0 ; j < n ; j++) {
					clone.rotate(*newTarget) ;
					clone.boost(tthrust) ;
					Collision* col = clone.collision(target) ;
					if(col != NULL) {
						target = newTarget ;
						thrust = tthrust ;
						break ;
					}
					else {
						clone.move(1.0) ;
					}
				}
			}
			delete col ;
		}

		Unit& CP = *target ;
		float angle = pod->diffAngle(CP) ;

		//we try to correct the angle
		//correction is calculated to counteract the angle between direction and speed
		float correction = 0. ;

		if(abs(angle) <= 18. && speed != 0 && thrust != 0) {
			correction = 2 * pod->diffSpeedAngle(CP) * thrust / speed ;
		}

		float returnAngle = pod->diffAngle(CP) + correction ;
		if(returnAngle > 18.0) returnAngle = 18.0 ;
		else if(returnAngle < -18.0) returnAngle = -18.0 ;

		return Move(returnAngle, thrust) ;
	}

	//block opponent
	Move computeBMove(Pod* pod, Pod* opp) {
		//pod is our blocker, opp is the opponent to block
		float dist2opp = pod->distanceSq(*opp) ;
		Point oppDest = Point(opp->getX() + opp->getVx(), opp->getY() + opp->getVy()) ;

		Point* target = NULL ;
		float thrust = 0 ;

		//on fonce sur l'adv à proximité
		if(dist2opp < 5000 * 5000) {
			target = &oppDest ;
			thrust = 200. ;
		}
		//si a proximité du CP de l'adv, on attend
		else if(pod->distanceSq(game->checkpoints[opp->getNextCheckpointId()]) < 2000 * 2000) {
			target = &oppDest ;
			thrust = 0. ;
		}
		else if(pod->distanceSq(game->checkpoints[opp->getNextCheckpointId() + 1 % game->getCheckpointCount()]) < 2000 * 2000) {
			target = &oppDest ;
			thrust = 0. ;
		}
		else if(pod->distanceSq(game->checkpoints[opp->getNextCheckpointId()]) < opp->distanceSq(game->checkpoints[opp->getNextCheckpointId()])){
			target = game->checkpoints + opp->getNextCheckpointId() ;
			thrust = 100 ;
		}
		else {
			target = game->checkpoints + (opp->getNextCheckpointId() + 1) % game->getCheckpointCount() ;
			thrust = 100 ;
		}

		return Move(pod->diffAngle(*target),thrust) ;
	}

	//TODO : à affiner
	bool checkCollision(Pod* p, const Move& m, Pod* o) {
		Pod pod = Pod(*p) ;
		pod.rotate(m.getAngle()) ;
		pod.boost(m.getThrust()) ;

		Pod opp = Pod(*o) ;
		opp.rotate(game->checkpoints[opp.getNextCheckpointId()]) ;
		opp.boost(200.);

		Collision* col = pod.collision(&opp) ;
		if(col != NULL) {
			delete col ;
			return true ;
		}
		return false ;
	}

public:
	Move* computeMoves() {
		Move* moves = new Move[2] ;

		//best opponent
		Pod* opp = std::max(oppPods,oppPods+1,
				[this](Pod* a,Pod* b) {
					return a->score(game->checkpoints + a->getNextCheckpointId()) < b->score(game->checkpoints + b->getNextCheckpointId()) ;
				}
		);

		bool firstTurn = myPods[0].getChecked() < game->getCheckpointCount() &&
				myPods[1].getChecked() < game->getCheckpointCount() ;

		float score0 = myPods[0].score(game->checkpoints + myPods[0].getNextCheckpointId()) ;
		float score1 = myPods[1].score(game->checkpoints + myPods[1].getNextCheckpointId()) ;

		if(firstTurn) {
			moves[0] = computeAMove(myPods) ;
			moves[1] = computeAMove(myPods + 1) ;
		}
		else if(score0 > score1) {
		//if(score0 > score1) {
			moves[0] = computeAMove(myPods) ;
			moves[1] = computeBMove(myPods + 1, opp) ;
		}
		else {
			moves[0] = computeBMove(myPods, opp) ;
			moves[1] = computeAMove(myPods + 1) ;
		}

		if(checkCollision(myPods, moves[0], oppPods) || checkCollision(myPods, moves[0], oppPods +1)) moves[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, moves[1], oppPods) || checkCollision(myPods+1, moves[1], oppPods+1)) moves[1].setThrust(-1.) ;

		return moves ;
	}

	Solution* computeSolution() {
		Solution* s = new Solution() ;

		//on crée une simulation, et on sauvegarde l'état courant du jeu
		Game simulation = Game(*game) ;
		Game* saved = this->game ;
		this->game = &simulation ;

		//simulation des tours consécutifs pour construction de la solution
		for(int i = 0 ; i < LENGTH ; i++) {
			std::unique_ptr<Move[]> m(this->computeMoves()) ;
			s->solution[i] = m[0] ;
			s->solution[i+LENGTH] = m[1] ;
			game->simulateNextTurn(m[0], m[1], Move(0,0), Move(0,0)) ;
		}

		//on restaure l'état initial de game
		this->game = saved ;

		return s;
	}
protected:
	Pod *myPods ;
	Pod *oppPods ;
};

//IA : the two pods try to race as fast as possible
class RaceIA : public SimpleIA {
public:
	RaceIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp=false) : SimpleIA(game, max_time, isOpp) {

	}
	~RaceIA() { }

	Move* computeMoves() {
		Move* mtab = new Move[2] ;
		mtab[0] = computeAMove(myPods) ;
		mtab[1] = computeAMove(myPods + 1) ;

		if(checkCollision(myPods, mtab[0], oppPods) || checkCollision(myPods, mtab[0], oppPods +1)) mtab[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, mtab[1], oppPods) || checkCollision(myPods+1, mtab[1], oppPods+1)) mtab[1].setThrust(-1.) ;

		return mtab ;
	}
};

class BlockIA : public SimpleIA {
public:
	BlockIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp=false) : SimpleIA(game, max_time, isOpp) {

	}
	~BlockIA() { }

	Move* computeMoves() {
		//best opponent
		Pod* opp = std::max(oppPods,oppPods+1,
				[this](Pod* a,Pod* b) {
					return a->score(game->checkpoints + a->getNextCheckpointId()) < b->score(game->checkpoints + b->getNextCheckpointId()) ;
				}
		);

		Move* mtab = new Move[2] ;
		mtab[0] = computeBMove(myPods, opp) ;
		mtab[1] = computeBMove(myPods + 1, opp) ;

		if(SimpleIA::checkCollision(myPods, mtab[0], oppPods) || checkCollision(myPods, mtab[0], oppPods +1)) mtab[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, mtab[1], oppPods) || checkCollision(myPods+1, mtab[1], oppPods+1)) mtab[1].setThrust(-1.) ;

		return mtab ;
	}
};

//IA : dummy IA that does nothing
class NullIA : public IA {
public:
	NullIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool playOpp=false) : IA(game) { }

	Solution* computeSolution() {
		return new Solution() ;
	}

	Move* computeMoves() {
		return new Move[2] ;
	}

};

class SAIA : public IA {
private:
	const int NUM_ITERATION = 100 ; //number of iteration at each T°
	const float alpha = 0.97 ; //temperature reduction at each iteration
	const int INITIAL_TEMP = 100 ; //initial T°
	const bool KEEP_BEST = true ;

	std::unique_ptr<Solution> bestSolution ;
	std::unique_ptr<Solution> currentSolution ;

	static bool acceptance(float oldValue, float newValue, float temp) {
		Random& r = Random::getInstance();
		float deltaE = newValue - oldValue ;

		if(deltaE >= 0) return true ;
		else {
			//boltzmann acceptance
			if(r.nextFloat() < exp(deltaE/temp)) return true ;
			else return false ;
		}
	}

public:
	SAIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp=false) : IA(game, max_time,isOpp) {
		bestSolution = nullptr ;
		currentSolution = nullptr ;
	}

	Solution* computeSolution() {
		float temp = INITIAL_TEMP ;

		//initial solution building & scoring
		if(!KEEP_BEST || bestSolution == nullptr) {
			//currentSolution.reset(Solution::randomSolution()) ;
			RaceIA ia = RaceIA(game) ;
			currentSolution.reset(ia.computeSolution()) ;

			this->scoreSolution<NullIA>(*currentSolution) ;
			bestSolution.reset(new Solution(*currentSolution)) ;
			std::cerr << "Initial solution at temp " << temp << " with score " << bestSolution->getScore() << std::endl;
		}
		else {
			bestSolution->shiftLeft() ;
			this->scoreSolution<NullIA>(*bestSolution) ;
			currentSolution.reset(new Solution(*bestSolution)) ;
			std::cerr << "Reusing best solution at temp " << temp << " with score " << bestSolution->getScore() << std::endl;
		}

		while(!timeOut()) {
			for(int i = 0 ; i < NUM_ITERATION && !timeOut() ; i++) {
				//creation of new solution
				std::unique_ptr<Solution> child(new Solution(*currentSolution)) ;
				child->mutate() ;

				//score the new solution against SimpleIA
				this->scoreSolution<NullIA>(*child);

				//is new solution accepted ?
				if(acceptance(currentSolution->getScore(), child->getScore(), temp)) {
					if(child->getScore() > bestSolution->getScore()) {
						bestSolution.reset(new Solution(*child)) ;
					}
					//switch to new accepted solution
					currentSolution.reset(child.release()) ;
				}
			}
			temp = temp * alpha ;
		}
		std::cerr << "Best solution at temp " << temp << " with score " << bestSolution->getScore() << std::endl;
		return bestSolution.release() ;
	}

	Move* computeMoves() {
		std::unique_ptr<Solution> s(this->computeSolution()) ;
		Move* moves = new Move[2] ;
		moves[0] = s->solution[0] ;
		moves[1] = s->solution[LENGTH] ;
		return moves ;
	}

};

int main()
{
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();

    Game game(laps,checkpointCount) ;
    //SimpleIA ia(&game) ;
    SAIA ia(&game) ;
    //SimpleIA ia(&game) ;

    // game loop
    while (1) {
    	ia.resetTimer();
    	game.readGame();

    	std::unique_ptr<Move[]> moves(ia.computeMoves()) ;

    	std::cerr << "Time : " << ia.getElapsedTime() << std::endl ;
 //   	std::cerr << "First move : " << currentSolution->solution[0] << std::endl ;
 ///   	std::cerr << "Second move : " << currentSolution->solution[LENGTH] << std::endl ;

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        // You have to output the target position
        // followed by the power (0 <= power <= 200)
        // i.e.: "x y power"

    	game.pods[0].Output(moves[0]);
    	game.pods[1].Output(moves[1]);
    }
}
