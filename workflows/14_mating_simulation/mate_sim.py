import random
import matplotlib.pyplot as plt

# -------------------------------
# Simulation Parameters
# -------------------------------
num_males = 2000
num_females = 2000

dt = 0.1        # time step

# Encounter and death rate parameters
alpha = 0.00001           # base encounter rate constant (per female per time unit)
death_rate_moving = 0.02   # death rate for a male while moving (per time unit)
death_rate_mating = 0.1    # chance a male dies during a mating event
death_rate_female = 0.001  # baseline death rate for females (per time unit)

# -------------------------------
# Initialize spiders
# -------------------------------
# Each individual is represented as a dictionary.
# We record an 'id', whether they are 'alive', their 'first_mating_time', and their 'sex'.
males = [{'id': i, 'alive': True, 'first_mating_time': None, 'sex': 'male'} for i in range(num_males)]
females = [{'id': i, 'alive': True, 'first_mating_time': None, 'sex': 'female'} for i in range(num_females)]

# -------------------------------
# Run the Simulation
# -------------------------------
# The simulation runs until every individual has either mated (first_mating_time recorded)
# or is dead.
time = 0.0
def simulation_incomplete():
    # Returns True if there is any individual that is still alive and has not mated.
    for spider in males + females:
        if spider['alive'] and spider['first_mating_time'] is None:
            return True
    return False

while simulation_incomplete():
    # ---------------------------
    # Update females (stationary)
    # ---------------------------
    for female in females:
        if female['alive']:
            if random.random() < death_rate_female * dt:
                female['alive'] = False

    # ---------------------------
    # Update males (moving and mating)
    # ---------------------------
    for male in males:
        if not male['alive']:
            continue  # skip dead males
        
        # Male might die while moving.
        if random.random() < death_rate_moving * dt:
            male['alive'] = False
            continue  # if he dies, he won’t mate in this time step
        
        # Determine the current list of available (alive) females.
        living_females = [f for f in females if f['alive']]
        
        # Only attempt a mate search if there is at least one female alive.
        if living_females:
            # The probability to encounter a female in this time step is:
            #   encounter_probability = (alpha * number_of_living_females) * dt
            encounter_probability = alpha * len(living_females) * dt
            if random.random() < encounter_probability:
                # Male has found a mate!
                # Record his first mating time if not already recorded.
                if male['first_mating_time'] is None:
                    male['first_mating_time'] = time

                # Randomly choose one of the living females.
                chosen_female = random.choice(living_females)
                # Record the female’s first mating time if this is her first mate.
                if chosen_female['first_mating_time'] is None:
                    chosen_female['first_mating_time'] = time

                # Now simulate the cost of mating for the male.
                # He might die as a result of the mating event.
                if random.random() < death_rate_mating:
                    male['alive'] = False

    # Increment simulation time.
    time += dt

# -------------------------------
# Summarize the Results and Write Report Files
# -------------------------------
# Collect first mating times for individuals that mated.
male_mating_times = [male['first_mating_time'] for male in males if male['first_mating_time'] is not None]
female_mating_times = [female['first_mating_time'] for female in females if female['first_mating_time'] is not None]

print("Male first mating times (for those that mated):")
print(male_mating_times)
print("\nFemale first mating times (for those that mated):")
print(female_mating_times)

# Write simulation parameters to a text file.
with open("simulation_parameters.txt", "w") as f:
    f.write("Simulation Parameters:\n")
    f.write(f"num_males: {num_males}\n")
    f.write(f"num_females: {num_females}\n")
    f.write(f"dt: {dt}\n")
    f.write(f"alpha: {alpha}\n")
    f.write(f"death_rate_moving: {death_rate_moving}\n")
    f.write(f"death_rate_mating: {death_rate_mating}\n")
    f.write(f"death_rate_female: {death_rate_female}\n")

# Write a tab-separated file with the first mating record for each successfully mated individual.
with open("first_mating_record.txt", "w") as f:
    # Header line
    f.write("ind_id\tfirst_mating_time\tsex\n")
    # Record for males
    for male in males:
        if male['first_mating_time'] is not None:
            f.write(f"male_{male['id']}\t{male['first_mating_time']}\t{male['sex']}\n")
    # Record for females
    for female in females:
        if female['first_mating_time'] is not None:
            f.write(f"female_{female['id']}\t{female['first_mating_time']}\t{female['sex']}\n")

# -------------------------------
# Plotting the Distribution of First Mating Times
# -------------------------------
plt.figure(figsize=(10, 5))
plt.hist(male_mating_times, bins=20, alpha=0.5, label="Males")
plt.hist(female_mating_times, bins=20, alpha=0.5, label="Females")
plt.xlabel("Time of first mating")
plt.ylabel("Count")
plt.title("Distribution of First Mating Times")
plt.legend()
plt.show()

