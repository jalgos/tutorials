/*The mongodb team put together a really nice tutorial here:
http://try.mongodb.org/
Other important resources:
-javascript:   http://www.w3schools.com/js/
-json:       http://www.w3schools.com/json/

//The tutorial doesn't cover any of the admin problematics. It's intended for an end user we need to query and fill up a database in MongoDB.
//Check the mongoadmin tutorial for information on that.

********************************************************************************
0)Introduction
The legend has that MongoDB was created to exactly mimic Google's BigTable's behaviour.
It was designed to overcome SQL's main drawbacks that hindered its scalability.

-SQL's query language is not a Turing complete language. Exotic queries are hard to achieve. 
-Schemas limit flexibility. Insert have to fit the schema and in case you need to insert more information on subsequent you need to change the schema
-Column based tables make insert slow

MongoDB differs from SQL on many features:
-It fully relies on javascript, so queries can take all forms
-It has no schemas. Tables are simply collections of json documents. Records in the same collection can totally differ from one another.
-Row based database, makes insert fast but queries slower.
-MongoDB answers the problem of scalibility with replication sets and sharding (dividing the database among many servers to distribute queries).
-MongoDB loose schemas don't allow for indexing and joining. To do such operations you need to do it manually using the language

Be aware that no tool is inherently better than the other (although the fact that SQL doesn't have a real query language is a huge limitation*). Everything depends on the kind of problem and data you have. It also depends on what constraints you have. Are you trying to be very flexible, do you need fast access to the db, or fast writing. If joining and indexing is a key feature you want to have in your database MongoDB is definitely not a good solution

*There are other great db solutions that are column based, allow replication and sharding and have a full query language.
 Anyway knowing MongoDB is definitely a must if you want to understand better what is a database and what purposes it fullfills. (It should be compulsory for any data scientist! ^^ )

********************************************************************************
1)Installation

First step is to install the mongo client:
easy install on debian like distros:
>sudo apt-get install mongodb


********************************************************************************
2)The mongo shell

let's now open a mongo shell:

>mongo

Will connect to the localserver. It is of little interest if you have no local server running or no database on the server.

However, one of mongo's really nice feature is that it can execute code javascript code.
The mongo shell can be use as a command by command javascript interpreter:*/

print("Hallo Welt");
2+1;
var a = "Hahahaha!";
a;

/*If you don't know javascript I strongly suggest that you do the w3school tutorial on js: http://www.w3schools.com/js/
Knowing js is not mandatory if you want to use mongo but it's a key when trying to do finer things like running exotic queries or using map/reduce.

********************************************************************************
3)Connecting to a database
If your mongodb admin gave you some logging credential,
you connect to the server hostname and the database test:

mongo -u username -p password -host hostname test
Depending on your user privileges you can run the following commands to explore the db:*/

show dbs;
/*It will list all the databases served as long with an estimated size in GB of the db. The command won't run for users that don't have the readAnyDatabase role (so in general it shouldn't work if you were given access to only one db)*/

show collections;
/*It will list all the collections within the database. Collections are the rough equivalent of tables in SQL. Since mongo and SQL are conceptually different we prefer using different words. */

/********************************************************************************
4)Manually create databases and collections
If you have the right privilege creating a database and a collection is straightforward.
(Advanced stuff)
If you want to know more about privileges: http://docs.mongodb.org/manual/reference/built-in-roles/
The privilege requires to create new databases or delete a database is clusterAdmin.
The privilege require to create or drop a collection is readWrite */
//Note that you can switch to a database that doesn't even exist, mongo will not complain about that
use adbthatdoesntexist;
//Mongo won't create anything until a collection is added
//let's go back to test where you should have the right permissions
use test;
//Requesting an undefined collection will have the same effect
db.collectionthatdoesntexist.find();

//Creating a collection is just done by inserting a new document
//use the insert command
//Here we insert a document in the mytest collection. (if mytest is already taken use a different name that is free)
db.mytest.insert({name:"Vincenzo Nibali",nickname:"Le requin de Messine", victoires:4});


//insert feeds on JSON documents as almost all the mongo functions do
//Now let's look at our collection
db.mytest.find();

//You noticed that mongo added a field "_id" in our document. It should look like this: "_id" : ObjectId("53d126e18662a20aa550ad2a")
//This is used by mongo to index records and ensure unicity. If I repeat the insert above I will have to records with two different "_id"
//A way to ensure unicity of records is to the "_id" yourself. 
db.mytest.insert({name:"Rafal Majka",skills:["climbing","tactics"], victoires:2,"_id":"majka"});


//Note that the document has a different form than the previous records we inserted. 
//Any JSON will do, there is no need for consistency among records. even nested JSONs work, or JSONs containing an array (as above)
//If you try to repeat the majka insert you'll get an error
db.mytest.insert({name:"Rafal Majka",skills:["climbing","tactics"], victoires:2,"_id":"majka"});


//Now we can write queries on our database. We simply give an argument to find. It must be a JSON document that partially match the records you try to find
db.mytest.find({name:"Vincenzo Nibali"});
db.mytest.find({skills:["climbing","tactics"]});
db.mytest.find({name:"Bernard Thévenet"});//empty result

//Now let's say we want to update our database. We use the update command.
//update takes two arguments, a query and a JSON document describing how to update the records
//If the JSON is a regular document it will just replace the records by this JSON
db.mytest.update({name:"Rafal Majka"},{name:"Rafal Majka",skills:["climbing","tactics"], victoires:2,"_id":"majka",Country:"Poland"});
//In that particular case it was a bit tedious because we just wanted to update one record and had to copy paste everything
//MongoDB has a set of operator to conveniently update records: http://docs.mongodb.org/manual/reference/operator/update/
//Let's look at $set
db.mytest.update({name:"Rafal Majka"},{$set:{favourite_food:"pasta"}});//this will set a new field in the record without deleting the previous record

//Now run this query:
db.mytest.update({name:"Vincenzo Nibali"},{$set:{Country:"Italy"}});
//And look at the result. Only one of our two Nibalis records has been updated.
//update's behaviour can be further modified in its third argument. Here you can set a bunch of flags. One of them is "multi". If set to true, the update will be applied to all the matching records.
db.mytest.update({name:"Vincenzo Nibali"},{$set:{Country:"Italy"}},{multi:true});

//Another flag is the upsert flag that will add a record if there is no match
db.mytest.update({name:"Bernard Thévenet"},{name:"Bernard Thévenet",favourite_food:"andouillettes"});//Nothing changed in the db

db.mytest.update({name:"Bernard Thévenet"},{name:"Bernard Thévenet",favourite_food:"andouillettes"},{upsert:true});//A new record in the db


//The function save can save us some time when updating. It behaves like insert if given a document with no id
db.mytest.save({name:"Bernard Hinault",nickname:"The Badger"})
//It will overwrite any existing record matching the _id otherwise. It behave like update(...,{upsert:true}) but there is no need to write a query
db.mytest.save({name:"Bernard Hinault",nickname:"The Badger",_id:"majka"})

//Deleting records
//records are deleted using the "remove" function. It takes a JSON document as argument. This JSON document is a query in the collection
db.mytest.remove({name:"Bernard Hinault"})


//Databases creations, collections creations, filling the databases won't be done through the mongo shell. There are a lot of API's in various languages (Python, Perl, C++, Java,...) and they all use the same commands (roughly)

//Now that we went through the basic commands we'll look into a "real" database to do more advanced queries
//Let's just clean up our mess:
db.mytest.drop()

/********************************************************************************
5)Database operations
//From now on we assume you are connected to the sporcast test database. Credential should have been given to you*/

//let's look at the content of the play by play database:
db.pbp.find();//Will create a "cursor" (http://docs.mongodb.org/manual/core/cursors/) to the pbp data. It will display the first N record and you can view all of them if you keep typing:
it;
//A cursor is a pointer to a record in the database. It is used to loop through the records
//getting the number of records is pretty easy
db.pbp.count();

//It is of course possible to make finer queries than just matching records exactly. There are operators and you can use regex. The exhaustive list can be found there: http://docs.mongodb.org/manual/reference/operator/query/

//Greater than:
db.pbp.find({date:{$gt:ISODate("2013-02-27T00:00:00Z")}}).count()//All records which date is greater than 20130227

//Regex
db.pbp.find({play:/lebron\s*james/i})//All records where LeBron James was involved (case insensitive search)

//queries can be combined easily
db.pbp.find({date:{$lt:ISODate("2010-02-27T00:00:00Z")},play:/lebron\s*james/i})

//you can sort results
db.pbp.find({date:{$lt:ISODate("2010-02-27T00:00:00Z")},play:/lebron\s*james/i}).sort({$play_id:-1})//results sorted by decreasing play_id

//If you make non trivial queries you can use the $where operator. $where's value is either a string evaluating to javascript code, or a function.

db.pbp.find({$where:"parseInt(this.period_number)*10==this.score_home-2"})//this is the collection record
db.pbp.find({$where:function(){var reg=/(\d+)\s*ft\s*jumper/i;
			       var ms = reg.exec(this.play);
			       if(ms===null)return false;
			       return parseInt(ms[1])>10;}})//return jumpers from at least 10ft



//MongoDB provides a way to loop through all the records with the forEach function.
//It is applied to a cursor, it takes a function (javascript of course) as argument. The function takes the current record as argument. 
db.pbp.find().limit(50).forEach(function(obj){print(obj.play)})//this will print the first 50 plays in the db

//forEach can be used in conjunction with save to modify the db

db.pbp.find().limit(50).forEach(function(obj){obj.pos="Top 50";db.pbp.save(obj);})
//forEach can typically be used for joining


//MongoDB has many functions that can be use to manipulate the collections. The exhaustive list is there: http://docs.mongodb.org/manual/reference/method/js-collection/
//Some won't work on replica set or sharded databases.

//distinct
//An interesting function to explore the database is distinct. Even though collections in mongo are schemaless, it's possible to request the unique occurences of a field in a collection. This is similar to the distinct function in SQL
var locations=db.pbp.distinct("location");
locations;
//distinct returns an array
locations.length;

//Another key function is "aggregate", it is similar to the "group by" function in SQL. http://docs.mongodb.org/manual/reference/method/db.collection.aggregate/
//aggregate takes an array as its first input and options as JSON doc as its second input

db.pbp.aggregate([{$group:{_id:"$location",count:{$sum:1}}}])//The number of records by arena

//Starting with MongoDB 2.6 Results of aggregate can directly be put to a new collection using $out
//db.pbp.aggregate([{$group:{_id:"$location",count:{$sum:1}}},{$out:"db.location"}])

/*********************************************************************************/

//6)MongoDB utility scripts
//MongoDB also offers a full array of shell scripts to monitor, query or import data into the db

//The mongoexport tool.
//You can export a collection to a text containing the records as json documents.
//replace hostname/username/password with the right credentials
mongoexport --db test --host hostname --username username --password password --collection pbp --out pbp.json 
//There is also a feature that allows you to export the data to a csv file. Even though MongoDB is schemaless, this is doable by specifying the fields you want. Querying for a field that is not in the collection won't generate an error.
//The command is similar to above but you need to specify the fields(columns) and add the --csv flag
mongoexport --db test --host hostname --username username --password password --csv --collection pbp --out pbp.csv --fields season,date,date_time,game_id,location,game_note,period,time,score_home,score_away,off_team,def_team,play,play_id,game_time,game_date_time,period_number,period_type 
//conversely you can import data from flat files using mongoimport
//create a test.json file containing something like this
---
{location:"Kombank Arena",team:"Кошаркашки клуб Партизан"}
{location:"Olympic Indoor Hall",team:"Παναθηναϊκός"}
---

mongoimport--db test --host hostname --username username --password password --collection my_test test.json 

//Now with a csv file
/*test.csv*/
---
location,team
Antares,Le Mans
Nokia Arena,מ.כ. מכבי אלקטרה תל-אביב
---

mongoimport --db test --host hostname --username username --password  --collection location --file test.csv  --headerline --type csv

//Other utilities include mongodump and mongorestore. It's used to create a flat file backup of a db or create a db from such a dump. It's a useful tool to copy data across servers that are not in the same replica set. But this is beyond the scope of this tutorial, more on that on the mongoadmin tutorial.
//MongoDB also provides a lot of scripts to monitor performance or usage, or to deal with admin stuff.



/*********************************************************************************/
//7)Map/Reduce
//I thought I would do an entire part about map reduce, since this tool is really popular now.
//Map reduce is a standardize divide and conquer approach to querying the db. Map/Reduce is the core concept of Hadoop, however Hadoop provides a way more general usage than MongoDB's.
//Resource:
//http://en.wikipedia.org/wiki/MapReduce
//http://docs.mongodb.org/manual/core/map-reduce/

//Map/Reduce is a standardized implementation of divide and conquer. 
//Map/Reduce has 2 main parts. Map and Reduce as you would have guessed.
//Map will produced a pair of a key and a value to be emitted to the reduce function.
//The value is what we want to compute in our data. In the particular case of MongoDB this value has to be a scalar (i.e. not an object or an array).
//The reduce function will take the keys and values emitted in a subdivision of the db. The algorithm is designed so values are grouped by the same key. So the keys vector argument to the reduce function will be a repetition of the same key.
//The reduce function is used to combine values together to return one scalar. 
//The result of reduce will then be sent to another reduce function until all the databased is covered.
//You can see the graph generated by Map/Reduce calls as a tree



//Example:
//We want to know the distinct keys (columns) of a collection
//The map function will get the keys of the current object and emit the keys as a comma separated string. Since we don't want to group by any key at first, the keys emitted will be "1"
//Remmeber the value emitted has to be a scalar.
var all_keys_map = function(){
    var keys = Object.keys(this);
    emit("1",keys.join(","))
}


//Function that add values of an array into a set (here a js object)
//We'll be used to get the unique keys
var merge_arrays = function(sett,arr){
    for(var i=0; i<arr.length; ++i) 
    {
        sett[arr[i]]=1;
    }
}

//The reduce function will transform the values into arrays and then it will get the unique occurences. Then it repackages the unique keys into a comma separated string
var all_keys_reduce = function(keys,values){
    var sett = {};
    for(var i=0;i<values.length;i++){
        var arr = values[i].split(",");
        merge_arrays(sett,arr);
    }
    var res = Object.keys(sett);
    return res.join(",");
}

//You can give map/reduce a finalize function that will be called after the reducing is finished
//In our case we want to transform the reduced value back into an array
var all_keys_finalize = function(keys,reducedValue){
    return reducedValue.split(",");
}

//This is how you call map/reduce
db.pbp.mapReduce(all_keys_map,all_keys_reduce,{out:{inline:1},scope:{merge_arrays:merge_arrays},finalize:all_keys_finalize})
//The first two arguments are the map and reduce function. The third argument is a JSON to parametrize mapReduces's behavior. 
//The out field gives the collection to which mongodb should write the result. If it is set to inline, MongoDB will return a js object to the shell. 
//The finalize function is specified in the option.
//Map and Reduce can't take additional arguments. To further parametrize your function you can use global parameters. By default map and reduce functions work in an empty global environment. The scope field is used to set the environment variable.
//In our example the merge_arrays function had to be set.

//To familiarize yourself with map/reduce look at this example and complete it by writing a function that can compute the mean and variance of any "column" by any group
//Moments:
var moment_map = function(){
    if(this[moment_key]!==null)	emit(this[moment_key],Math.pow(this[mfield],morder));
}

var moment_reduce = function(keys,values){
    var res = 0;
    for(var i=0;i<values.length;i++){
	res+=values[i];
    }
    return res;
}

db.pbp.mapReduce(moment_map,moment_reduce,{out:{inline:1},scope:{moment_key:"def_team",mfield:"score_home",morder:2}})//sum of score_home^2 by def_team

//Map Reduce is particularly efficient on sharded database. Each reduce fonction will concurrently run on the different servers, hence multiplying the speed of execution by the number of shards in your db.

/*********************************************************************/

//This tutorial just tackles the basics of mongodb. If you want to go further just visit the mongodb site. Remember also that Google/stackoverflow is always your friend.
