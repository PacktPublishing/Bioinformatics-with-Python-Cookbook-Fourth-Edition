#!/usr/bin/env python3
"""
Complete deployable bioinformatics workflow using AWS Step Functions
This creates a simple workflow with one mocked FASTQ quality analysis task
"""

# Import Libraries
import boto3
import json
import time
import zipfile
import io
from botocore.exceptions import ClientError

# Configuration
REGION = 'us-west-2'  # Change to your preferred region
LAMBDA_FUNCTION_NAME = 'fastq-quality-analysis'
STATE_MACHINE_NAME = 'simple-bioinformatics-pipeline'
IAM_ROLE_NAME = 'BioinformaticsWorkflowRole'

def create_lambda_function_code():
    """Create the Lambda function code for FASTQ quality analysis"""
    lambda_code = '''
import json
import random
import time

def lambda_handler(event, context):
    """
    Mock FASTQ quality analysis function
    Simulates processing time and returns mock results
    """
    
    # Extract input parameters
    sample_id = event.get('sampleId', 'unknown')
    fastq_file = event.get('fastqFile', '')
    
    print(f"Starting quality analysis for sample: {sample_id}")
    print(f"Processing file: {fastq_file}")
    
    # Simulate processing time (2-5 seconds)
    processing_time = random.uniform(2, 5)
    time.sleep(processing_time)
    
    # Generate mock quality metrics
    total_reads = random.randint(1000000, 5000000)
    quality_score = round(random.uniform(25, 40), 2)
    gc_content = round(random.uniform(40, 60), 2)
    duplicate_rate = round(random.uniform(5, 15), 2)
    
    # Determine quality status
    if quality_score >= 30 and duplicate_rate <= 10:
        status = "PASS"
    elif quality_score >= 25 and duplicate_rate <= 15:
        status = "WARNING" 
    else:
        status = "FAIL"
    
    # Create results
    results = {
        "sampleId": sample_id,
        "inputFile": fastq_file,
        "status": status,
        "metrics": {
            "totalReads": total_reads,
            "averageQualityScore": quality_score,
            "gcContent": gc_content,
            "duplicateRate": duplicate_rate,
            "processingTimeSeconds": round(processing_time, 2)
        },
        "recommendations": []
    }
    
    # Add recommendations based on results
    if quality_score < 30:
        results["recommendations"].append("Consider quality trimming")
    if duplicate_rate > 10:
        results["recommendations"].append("High duplicate rate detected")
    if gc_content < 45 or gc_content > 55:
        results["recommendations"].append("Unusual GC content detected")
    
    print(f"Analysis complete for {sample_id}: {status}")
    
    return {
        'statusCode': 200,
        'body': results
    }
'''
    
    # Create ZIP file for Lambda deployment
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        zip_file.writestr('lambda_function.py', lambda_code)
    
    return zip_buffer.getvalue()

def create_state_machine_definition():
    """Create the Step Functions state machine definition"""
    return {
        "Comment": "Simple Bioinformatics Pipeline - FASTQ Quality Analysis",
        "StartAt": "ValidateInput",
        "States": {
            "ValidateInput": {
                "Type": "Pass",
                "Comment": "Validate input parameters",
                "Parameters": {
                    "sampleId.$": "$.sampleId",
                    "fastqFile.$": "$.fastqFile",
                    "timestamp.$": "$$.Execution.StartTime"
                },
                "Next": "FastqQualityAnalysis"
            },
            "FastqQualityAnalysis": {
                "Type": "Task",
                "Resource": "arn:aws:states:::lambda:invoke",
                "Parameters": {
                    "FunctionName": LAMBDA_FUNCTION_NAME,
                    "Payload": {
                        "sampleId.$": "$.sampleId",
                        "fastqFile.$": "$.fastqFile"
                    }
                },
                "ResultPath": "$.qualityAnalysis",
                "Next": "CheckQualityResults",
                "Retry": [
                    {
                        "ErrorEquals": ["States.TaskFailed"],
                        "IntervalSeconds": 2,
                        "MaxAttempts": 2,
                        "BackoffRate": 2.0
                    }
                ],
                "Catch": [
                    {
                        "ErrorEquals": ["States.ALL"],
                        "Next": "AnalysisFailure",
                        "ResultPath": "$.error"
                    }
                ]
            },
            "CheckQualityResults": {
                "Type": "Choice",
                "Choices": [
                    {
                        "Variable": "$.qualityAnalysis.Payload.body.status",
                        "StringEquals": "PASS",
                        "Next": "AnalysisSuccess"
                    },
                    {
                        "Variable": "$.qualityAnalysis.Payload.body.status", 
                        "StringEquals": "WARNING",
                        "Next": "AnalysisWarning"
                    }
                ],
                "Default": "AnalysisFailure"
            },
            "AnalysisSuccess": {
                "Type": "Pass",
                "Parameters": {
                    "result": "SUCCESS",
                    "message": "FASTQ quality analysis passed all checks",
                    "sampleId.$": "$.sampleId",
                    "qualityMetrics.$": "$.qualityAnalysis.Payload.body.metrics"
                },
                "End": True
            },
            "AnalysisWarning": {
                "Type": "Pass", 
                "Parameters": {
                    "result": "WARNING",
                    "message": "FASTQ quality analysis completed with warnings",
                    "sampleId.$": "$.sampleId",
                    "qualityMetrics.$": "$.qualityAnalysis.Payload.body.metrics",
                    "recommendations.$": "$.qualityAnalysis.Payload.body.recommendations"
                },
                "End": True
            },
            "AnalysisFailure": {
                "Type": "Pass",
                "Parameters": {
                    "result": "FAILURE", 
                    "message": "FASTQ quality analysis failed",
                    "sampleId.$": "$.sampleId",
                    "error.$": "$.error"
                },
                "End": True
            }
        }
    }

def create_iam_role():
    """Create IAM role for Step Functions and Lambda"""
    iam = boto3.client('iam', region_name=REGION)
    
    # Trust policy for the role
    trust_policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {
                    "Service": [
                        "states.amazonaws.com",
                        "lambda.amazonaws.com"
                    ]
                },
                "Action": "sts:AssumeRole"
            }
        ]
    }
    
    # Permission policy
    permission_policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": [
                    "lambda:InvokeFunction",
                    "logs:CreateLogGroup",
                    "logs:CreateLogStream", 
                    "logs:PutLogEvents"
                ],
                "Resource": "*"
            }
        ]
    }
    
    try:
        # Create role
        role_response = iam.create_role(
            RoleName=IAM_ROLE_NAME,
            AssumeRolePolicyDocument=json.dumps(trust_policy),
            Description='Role for bioinformatics Step Functions workflow'
        )
        
        # Attach basic Lambda execution policy
        iam.attach_role_policy(
            RoleName=IAM_ROLE_NAME,
            PolicyArn='arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole'
        )
        
        # Create and attach custom policy
        policy_response = iam.create_policy(
            PolicyName=f'{IAM_ROLE_NAME}Policy',
            PolicyDocument=json.dumps(permission_policy),
            Description='Custom policy for bioinformatics workflow'
        )
        
        iam.attach_role_policy(
            RoleName=IAM_ROLE_NAME,
            PolicyArn=policy_response['Policy']['Arn']
        )
        
        print(f"✅ Created IAM role: {IAM_ROLE_NAME}")
        return role_response['Role']['Arn']
        
    except ClientError as e:
        if e.response['Error']['Code'] == 'EntityAlreadyExists':
            role = iam.get_role(RoleName=IAM_ROLE_NAME)
            print(f"✅ IAM role already exists: {IAM_ROLE_NAME}")
            return role['Role']['Arn']
        else:
            raise

def create_lambda_function(role_arn):
    """Create the Lambda function"""
    lambda_client = boto3.client('lambda', region_name=REGION)
    
    try:
        response = lambda_client.create_function(
            FunctionName=LAMBDA_FUNCTION_NAME,
            Runtime='python3.9',
            Role=role_arn,
            Handler='lambda_function.lambda_handler',
            Code={'ZipFile': create_lambda_function_code()},
            Description='Mock FASTQ quality analysis for bioinformatics pipeline',
            Timeout=30,
            MemorySize=128
        )
        print(f"✅ Created Lambda function: {LAMBDA_FUNCTION_NAME}")
        return response['FunctionArn']
        
    except ClientError as e:
        if e.response['Error']['Code'] == 'ResourceConflictException':
            # Update existing function
            lambda_client.update_function_code(
                FunctionName=LAMBDA_FUNCTION_NAME,
                ZipFile=create_lambda_function_code()
            )
            response = lambda_client.get_function(FunctionName=LAMBDA_FUNCTION_NAME)
            print(f"✅ Updated existing Lambda function: {LAMBDA_FUNCTION_NAME}")
            return response['Configuration']['FunctionArn']
        else:
            raise

def create_step_function(role_arn):
    """Create the Step Functions state machine"""
    stepfunctions = boto3.client('stepfunctions', region_name=REGION)
    
    definition = create_state_machine_definition()
    
    try:
        response = stepfunctions.create_state_machine(
            name=STATE_MACHINE_NAME,
            definition=json.dumps(definition),
            roleArn=role_arn,
            type='STANDARD'
        )
        print(f"✅ Created State Machine: {STATE_MACHINE_NAME}")
        return response['stateMachineArn']
        
    except ClientError as e:
        if e.response['Error']['Code'] == 'StateMachineAlreadyExists':
            # Update existing state machine
            machines = stepfunctions.list_state_machines()
            machine_arn = None
            for machine in machines['stateMachines']:
                if machine['name'] == STATE_MACHINE_NAME:
                    machine_arn = machine['stateMachineArn']
                    break
            
            if machine_arn:
                stepfunctions.update_state_machine(
                    stateMachineArn=machine_arn,
                    definition=json.dumps(definition),
                    roleArn=role_arn
                )
                print(f"✅ Updated existing State Machine: {STATE_MACHINE_NAME}")
                return machine_arn
        else:
            raise

def run_workflow_example(state_machine_arn):
    """Run an example workflow execution"""
    stepfunctions = boto3.client('stepfunctions', region_name=REGION)
    
    # Example input data
    input_data = {
        "sampleId": "SAMPLE_001",
        "fastqFile": "s3://my-bioinformatics-bucket/samples/SAMPLE_001_R1.fastq.gz"
    }
    
    try:
        response = stepfunctions.start_execution(
            stateMachineArn=state_machine_arn,
            name=f"test-execution-{int(time.time())}",
            input=json.dumps(input_data)
        )
        
        execution_arn = response['executionArn']
        print(f"✅ Started workflow execution: {execution_arn}")
        
        # Wait for completion
        print("⏳ Waiting for execution to complete...")
        while True:
            status = stepfunctions.describe_execution(executionArn=execution_arn)
            if status['status'] in ['SUCCEEDED', 'FAILED', 'TIMED_OUT', 'ABORTED']:
                break
            time.sleep(2)
        
        print(f"✅ Execution completed with status: {status['status']}")
        
        if status['status'] == 'SUCCEEDED':
            output = json.loads(status['output'])
            print("📊 Execution Results:")
            print(json.dumps(output, indent=2))
        
        return status
        
    except Exception as e:
        print(f"❌ Error running workflow: {e}")
        return None

def main():
    """Main deployment function"""
    print("🚀 Deploying Simple Bioinformatics Workflow...")
    print(f"📍 Region: {REGION}")
    
    try:
        # Step 1: Create IAM role
        print("\n1️⃣ Creating IAM role...")
        role_arn = create_iam_role()
        
        # Wait for role propagation
        print("⏳ Waiting for IAM role propagation...")
        time.sleep(10)
        
        # Step 2: Create Lambda function
        print("\n2️⃣ Creating Lambda function...")
        lambda_arn = create_lambda_function(role_arn)
        
        # Step 3: Create Step Functions state machine
        print("\n3️⃣ Creating Step Functions state machine...")
        state_machine_arn = create_step_function(role_arn)
        
        # Step 4: Run example workflow
        print("\n4️⃣ Running example workflow...")
        run_workflow_example(state_machine_arn)
        
        print(f"\n🎉 Deployment completed successfully!")
        print(f"📋 State Machine ARN: {state_machine_arn}")
        print(f"🔧 Lambda Function ARN: {lambda_arn}")
        
        print(f"\n💡 To run the workflow manually:")
        print(f"aws stepfunctions start-execution \\")
        print(f"  --state-machine-arn {state_machine_arn} \\")
        print(f"  --input '{json.dumps({'sampleId': 'TEST_SAMPLE', 'fastqFile': 's3://bucket/file.fastq.gz'})}'")
        
    except Exception as e:
        print(f"❌ Deployment failed: {e}")
        raise

def cleanup():
    """Clean up all created resources"""
    print("🧹 Cleaning up resources...")
    
    # Delete state machine
    try:
        stepfunctions = boto3.client('stepfunctions', region_name=REGION)
        machines = stepfunctions.list_state_machines()
        for machine in machines['stateMachines']:
            if machine['name'] == STATE_MACHINE_NAME:
                stepfunctions.delete_state_machine(stateMachineArn=machine['stateMachineArn'])
                print(f"🗑️ Deleted state machine: {STATE_MACHINE_NAME}")
                break
    except Exception as e:
        print(f"⚠️ Error deleting state machine: {e}")
    
    # Delete Lambda function
    try:
        lambda_client = boto3.client('lambda', region_name=REGION)
        lambda_client.delete_function(FunctionName=LAMBDA_FUNCTION_NAME)
        print(f"🗑️ Deleted Lambda function: {LAMBDA_FUNCTION_NAME}")
    except Exception as e:
        print(f"⚠️ Error deleting Lambda function: {e}")
    
    # Delete IAM role and policies
    try:
        iam = boto3.client('iam', region_name=REGION)
        
        # Detach policies
        attached_policies = iam.list_attached_role_policies(RoleName=IAM_ROLE_NAME)
        for policy in attached_policies['AttachedPolicies']:
            iam.detach_role_policy(RoleName=IAM_ROLE_NAME, PolicyArn=policy['PolicyArn'])
            if 'BioinformaticsWorkflowRolePolicy' in policy['PolicyName']:
                iam.delete_policy(PolicyArn=policy['PolicyArn'])
        
        # Delete role
        iam.delete_role(RoleName=IAM_ROLE_NAME)
        print(f"🗑️ Deleted IAM role: {IAM_ROLE_NAME}")
    except Exception as e:
        print(f"⚠️ Error deleting IAM role: {e}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == "cleanup":
        cleanup()
    else:
        main()
